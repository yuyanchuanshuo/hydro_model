program main
    implicit none
    integer :: i, j, number, t_all, dt, t,j_out
    real, dimension(:, :), allocatable :: A, Q, Q0, A0, H, R
    real, dimension(:), allocatable :: n, dx, ZB
    integer :: file_unit  ! 添加文件单元号变量
    integer :: file_unit_all    ! 添加专门输出第2000时刻所有断面的文件

    ! 新增断面数据结构
    type :: section_data
        character(len=256) :: id
        real :: start_distance
        real, dimension(:), allocatable :: distances
        real, dimension(:), allocatable :: elevations
    endtype section_data

    type :: hydraulic_data
        character(len=256) :: section_id
        real, dimension(:, :), allocatable :: results  ! 水位、面积、宽度、湿周、水力半径
    endtype hydraulic_data

    type(section_data), dimension(:), allocatable :: sections
    type(hydraulic_data), dimension(:), allocatable :: section_hydraulics
    integer :: n_sections

    ! 读取断面数据
    call read_section_data('data_input/cross_section_data.csv', sections, n_sections)
    allocate(section_hydraulics(n_sections))

    ! 计算水力特性
    do i = 1, n_sections
        section_hydraulics(i)%section_id = sections(i)%id
        call calculate_hydraulic_properties(sections(i), section_hydraulics(i)%results)
    enddo

    ! 新增：输出到文件
    open(newunit=file_unit, file='data_output/section_hydraulic_properties.txt', status='replace', action='write')
    write(file_unit, '(A)') '断面ID 水位(m) 面积(m2) 水面宽度(m) 湿周(m) 水力半径(m)'

    do i = 1, n_sections
        do j = 1, size(section_hydraulics(i)%results, 1)
            write(file_unit, '(A10, 5F15.6)') trim(section_hydraulics(i)%section_id), &
                section_hydraulics(i)%results(j, 1), &  ! 水位
                section_hydraulics(i)%results(j, 2), &  ! 面积
                section_hydraulics(i)%results(j, 3), &  ! 水面宽度
                section_hydraulics(i)%results(j, 4), &  ! 湿周
                section_hydraulics(i)%results(j, 5)     ! 水力半径
        enddo
    enddo
    close(file_unit)

    number = n_sections !假设有100个断面，待修改
    t_all = 100 !计算总时段
    dt = 3 ! 计算间隔10s

    ! 分配数组内存
    allocate(A(number, t_all), A0(number, t_all), Q(number, t_all), Q0(number, t_all), R(number, t_all), H(number, t_all))!定义以断面为行，时间为列的矩阵
    allocate(n(number), dx(number), ZB(number))

    Q0(:, 1) = 100 !假设每一个断面的第一时刻流量均为100m3/s
    Q(:, 1) = 100 !假设每一个断面的第一时刻流量均为100m3/s
    A0(:, 1) = 200 !假设每一个断面的第一时刻面积均为1000m2
    A(:, 1) = 200 !假设每一个断面的第一时刻面积均为1000m2
    A(number, :) = 100 ! 假设最后一个断面的所有时刻过水面积都是200

    n = 0.035 ! 假设每一个断面的糙率均为0.035

    do i = 2, number
        if(i == 1) then
            dx(i) = 0
            ZB(i) = minval(sections(i)%elevations)
        else
            dx(i) = abs(sections(i)%start_distance - sections(i - 1)%start_distance)
            ZB(i) = minval(sections(i)%elevations)
        endif
    enddo

    do i = 1, t_all ! 定义上边界入流，这个以后从文件中读取
        Q0(1, i) = 100 + (150 - 100) * real(i - 1) / (t_all - 1)
        Q(1, i) = 100 + (150 - 100) * real(i - 1) / (t_all - 1)
    enddo

    j_out = 100
    do t = 1, t_all - 1

        call dis_cal(number, n, dx, dt, Q0(:, t), ZB, A(:, t), Q(:, t + 1), H(:, t))  ! 使用call调用子例程
        Q0(2:number, t + 1) = Q(2:number, t + 1) ! 定义这个断面输入的Q
        Q(1, t + 1) = Q0(1, t + 1)
        A(:, t + 1) = area_cal(number, dt, Q(:, t + 1), A0(:, t), dx)
        A0(:, t + 1) = A(:, t + 1) ! 定义这个断面输入的A

        print *, '第', t, '时刻', j_out, '个断面的流量为：', Q(j_out, t), A(j_out, t), H(j_out, t)

    enddo

    ! 关闭文件
    close(file_unit)
    close(file_unit_all)

    ! 在程序结束前添加内存释放语句
    deallocate(A, Q, Q0, A0, R, H, n, dx, ZB)

contains
    function area_cal(n_size, dt_in, Q_in, A_in, dx_in) result(A_out)
        implicit none
        integer :: ii
        integer, intent(in) :: n_size
        integer, intent(in) :: dt_in
        real, intent(in), dimension(n_size) :: Q_in, A_in, dx_in
        real, dimension(n_size) :: A_out

        Do ii = 1, n_size - 1
            A_out(ii) = A_in(ii) - 2 * dt_in * (Q_in(ii + 1) - Q_in(ii)) / ((dx_in(ii + 1) + dx_in(ii)) / 2)
        enddo
        ! print *, A_out
    endfunction area_cal

    ! 添加处理起点距的函数
    function parse_distance(distance_str) result(distance_val)
        character(len=*), intent(in) :: distance_str
        real :: distance_val
        integer :: plus_pos, part1, part2

        plus_pos = index(distance_str, '+')
        if(plus_pos > 0) then
            read(distance_str(1:plus_pos - 1), *) part1
            read(distance_str(plus_pos + 1:), *) part2
            distance_val = real(part1 * 1000 + part2)
        else
            read(distance_str, *) distance_val
        endif
    endfunction parse_distance

    subroutine dis_cal(n_size, n_in, dx_in, dt_in, Q_in, ZB_in, A_in, Q_out, h_out) !函数只能返回一个结果变量，改用子例程
        implicit none
        integer :: iii
        integer, intent(in) :: n_size, dt_in
        real, intent(in), dimension(n_size) :: dx_in, n_in
        real, dimension(n_size) :: uu, Q_in, H_in, ZB_in, A_in, R_in, B, wetted_perimeter,ZZ
        real, dimension(n_size) :: Q_out, h_out
        real :: XDX, HDX, FRI, alfa = 1 !alfa动量修正系数，一般取1

        Do iii = 1, n_size
            call find_hydraulic_properties(section_hydraulics(iii)%results, A_in(iii), ZZ(iii), B(iii), wetted_perimeter(iii), R_in(iii))
            H_in(iii) = ZZ(iii) - ZB_in(iii)
        enddo

        h_out = H_in

        uu(1) = Q_in(1) / A_in(1)

        Do iii = 2, n_size
            uu(iii) = Q_in(iii) / A_in(iii)
            if(uu(iii) > 0) then ! 计算对流项
                if(uu(iii - 1) > 0) then
                    XDX = alfa * (uu(iii) * Q_in(iii) - uu(iii - 1) * Q_in(iii - 1)) / ((dx_in(iii) + dx_in(iii - 1)) / 2)
                else
                    XDX = alfa * (uu(iii) * Q_in(iii) - 0) / ((dx_in(iii) + dx_in(iii - 1)) / 2)
                endif
            elseif(uu(iii) < 0) then
                if(uu(iii + 1) > 0) then
                    XDX = alfa * (0 - uu(iii) * Q_in(iii)) / ((dx_in(iii + 1) + dx_in(iii)) / 2)
                else
                    XDX = alfa * (uu(iii + 1) * Q_in(iii + 1) - uu(iii) * Q_in(iii)) / ((dx_in(iii + 1) + dx_in(iii)) / 2)
                endif
            endif

            HDX = 9.8 * (H_in(iii) + ZB_in(iii) - H_in(iii - 1) - ZB_in(iii - 1)) / dx_in(iii) * & ! 使用续行符&表示续行
                  (A_in(iii) + A_in(iii - 1)) / 2

            FRI = 9.8 * ABS(uu(iii)) * ((n_in(iii) + n_in(iii - 1)) / 2)**2 / & ! 计算FRI，除数的那一个
                  (2 * ((R_in(iii) + R_in(iii - 1)) / 2)**1.333)

            Q_out(iii) = (Q_in(iii) - 2.*dt_in * XDX - 2.*dt_in * HDX - 2.*dt_in * FRI * Q_in(iii)) / (1.+2.*dt_in * FRI)

        enddo

    endsubroutine dis_cal 

    ! 添加断面数据读取子程序
    subroutine read_section_data(filename, sections_out, n_sections_out)
        character(len=*), intent(in) :: filename
        type(section_data), dimension(:), allocatable, intent(out) :: sections_out
        integer, intent(out) :: n_sections_out
        integer :: ios, section_num, ii, jj, n_points
        integer :: total_lines = 0 ! 文件总行数
        character(len=256) :: line, start_distance_str
        integer, dimension(:), allocatable :: id_line_numbers  ! 新增：存储ID行行号

        ! 打开文件并读取数据
        open(newunit=file_unit, file=filename, status='old', action='read')
        ! 第一遍读取确定断面数量
        do
            read(file_unit, '(A)', iostat=ios) line
            if(ios /= 0) exit
            total_lines = total_lines + 1
        enddo
        rewind(file_unit)

        ! 分配ID行行号数组
        allocate(id_line_numbers(total_lines))  ! 最大可能数量
        n_sections_out = 0

        do ii = 1, total_lines
            read(file_unit, '(A)') line
            ! 判断是否是断面ID行(行首是字符串)
            if(verify(line(1:1), '0123456789') > 0) then  ! 简化判断方法
                n_sections_out = n_sections_out + 1
                id_line_numbers(n_sections_out) = ii
            endif
        enddo

        ! 压缩数组，移除未使用的0元素
        if(n_sections_out < total_lines) then
            id_line_numbers = pack(id_line_numbers, id_line_numbers /= 0)
        endif

        rewind(file_unit) ! 回到文件起点

        ! 分配断面数据内存
        allocate(sections_out(n_sections_out))
        section_num = 0
        jj = 0
        ! 第二遍读取：逐行读取，判断是ID行，n就+1
        do ii = 1, total_lines

            if(ii == id_line_numbers(section_num + 1)) then
                if(ii /= 1) then
                    ! 调整数组大小为实际读取的点数
                    sections_out(section_num)%distances = sections_out(section_num)%distances(1:n_points)
                    sections_out(section_num)%elevations = sections_out(section_num)%elevations(1:n_points)
                endif

                section_num = section_num + 1
                ! 解析断面ID和起点距
                read(file_unit, '(A)', iostat=ios) line  ! 读取整行到line变量
                read(line, *) sections_out(section_num)%id, start_distance_str! 按空格分隔读取多个变量
                sections_out(section_num)%start_distance = parse_distance(start_distance_str)
                n_points = 0
                jj = jj + 1
            else
                ! 确保在使用前分配了数组内存
                if(.not. allocated(sections_out(section_num)%distances)) then
                    allocate(sections_out(section_num)%distances(1000))
                    allocate(sections_out(section_num)%elevations(1000))
                endif
                n_points = n_points + 1
                read(file_unit, '(A)', iostat=ios) line  ! 读取整行到line变量
                read(line, *) sections_out(section_num)%distances(n_points), sections_out(section_num)%elevations(n_points)
                jj = jj + 1

                ! 最后一行读取后调整section的大小
                if(ii == total_lines) then
                    ! 调整数组大小为实际读取的点数
                    sections_out(section_num)%distances = sections_out(section_num)%distances(1:n_points)
                    sections_out(section_num)%elevations = sections_out(section_num)%elevations(1:n_points)
                endif
            endif

        enddo
        close(file_unit)

        ! 添加断面倒序判断逻辑，确保最后一个断面最低
        if(minval(sections_out(1)%elevations) < minval(sections_out(n_sections_out)%elevations)) then
            sections_out = sections_out(n_sections_out:1:-1)  ! 倒序断面数组
        endif

        ! 新增：输出断面详细数据到文件
        open(newunit=file_unit, file='data_output/section_details.txt', status='replace', action='write')
        write(file_unit, '(A)') '断面序号 断面ID 起点距(m) 点号 距离(m) 高程(m)'

        do ii = 1, n_sections_out
            do jj = 1, size(sections_out(ii)%distances)
                write(file_unit, '(I10, A10, F15.3, I10, 2F15.3)') ii, trim(sections_out(ii)%id), &
                    sections_out(ii)%start_distance, &
                    jj, &
                    sections_out(ii)%distances(jj), &
                    sections_out(ii)%elevations(jj)
            enddo
        enddo
        close(file_unit)

    endsubroutine read_section_data

    ! 添加水力特性计算子程序
    subroutine calculate_hydraulic_properties(section_in, results_matrix)
        type(section_data), intent(in) :: section_in
        real, dimension(:, :), allocatable, intent(out) :: results_matrix
        real :: min_elev, max_elev, current_h
        integer :: iii, jjj, n_levels, n_points_in
        real :: left_dist, right_dist, left_elev, right_elev
        real :: segment_width, segment_area, segment_perimeter
        real, dimension(:), allocatable :: areas, widths, wetted_perimeters

        n_points_in = size(section_in%distances)
        min_elev = minval(section_in%elevations)
        max_elev = maxval(section_in%elevations)

        ! 设置水位范围：从最低点到最高点+10m，步长0.1m
        n_levels = int((max_elev + 10.0 - min_elev) / 0.1) + 1
        allocate(results_matrix(n_levels, 5))  ! 4列：水位、面积、宽度、湿周、水力半径
        allocate(areas(n_levels))
        allocate(widths(n_levels))
        allocate(wetted_perimeters(n_levels))

        ! 初始化水位数组
        do iii = 1, n_levels
            results_matrix(iii, 1) = min_elev + (iii - 1) * 0.1
        enddo

        ! 计算每个水位对应的水力特性
        do iii = 1, n_levels
            current_h = results_matrix(iii, 1)
            areas(iii) = 0.0
            widths(iii) = 0.0
            wetted_perimeters(iii) = 0.0

            ! 遍历所有断面点
            do jjj = 2, n_points_in
                left_dist = section_in%distances(jjj - 1)
                right_dist = section_in%distances(jjj)
                left_elev = section_in%elevations(jjj - 1)
                right_elev = section_in%elevations(jjj)

                ! 计算当前水位下的有效宽度
                if(current_h > min(left_elev, right_elev)) then
                    segment_width = right_dist - left_dist

                    ! 计算梯形面积
                    if(current_h >= max(left_elev, right_elev)) then
                        segment_area = segment_width * (current_h - (left_elev + right_elev) / 2.0)
                        segment_perimeter = segment_width + &
                                            sqrt((right_dist - left_dist)**2 + (right_elev - left_elev)**2)
                    else
                        ! 部分淹没情况(三角形)
                        if(left_elev < right_elev) then
                            segment_width = segment_width * (current_h - left_elev) / (right_elev - left_elev)
                            segment_area = 0.5 * segment_width * (current_h - left_elev)
                            segment_perimeter = sqrt(segment_width**2 + (current_h - left_elev)**2)
                        else
                            segment_width = segment_width * (current_h - right_elev) / (left_elev - right_elev)
                            segment_area = 0.5 * segment_width * (current_h - right_elev)
                            segment_perimeter = sqrt(segment_width**2 + (current_h - right_elev)**2)
                        endif
                    endif
                    widths(iii) = widths(iii) + segment_width
                    areas(iii) = areas(iii) + segment_area
                    wetted_perimeters(iii) = wetted_perimeters(iii) + segment_perimeter
                endif
            enddo

            ! 填充结果矩阵
            results_matrix(iii, 2) = areas(iii)                 ! 面积
            results_matrix(iii, 3) = widths(iii)                ! 水面宽度
            results_matrix(iii, 4) = wetted_perimeters(iii)     ! 湿周

            ! 计算水力半径
            if(wetted_perimeters(iii) > 0.0) then
                results_matrix(iii, 5) = areas(iii) / wetted_perimeters(iii)
            else
                results_matrix(iii, 5) = 0.0
            endif
        enddo

        deallocate(areas, widths, wetted_perimeters)

    endsubroutine calculate_hydraulic_properties

    ! 修改后的查找函数，直接使用已计算的水力特性数据
    subroutine find_hydraulic_properties(result_in, A_in, h_out, B_out, wp_out, R_out)
        real, dimension(:, :), allocatable, intent(in) :: result_in
        real, intent(in) :: A_in
        real, intent(out) :: h_out, B_out, wp_out, R_out
        integer :: ii

        ! 直接使用已计算的水力特性结果
        do ii = 1, size(result_in, 1) - 1
            if(result_in(ii, 2) <= A_in .and. &
               result_in(ii + 1, 2) >= A_in) then
                ! 线性插值
                h_out = result_in(ii, 1) + (A_in - result_in(ii, 2)) * (result_in(ii + 1, 1) - result_in(ii, 1)) / &
                              (result_in(ii + 1, 2) - result_in(ii, 2))
                B_out = result_in(ii, 3) + (A_in - result_in(ii, 2)) * (result_in(ii + 1, 3) - result_in(ii, 3)) / &
                        (result_in(ii + 1, 2) - result_in(ii, 2))
                wp_out = result_in(ii, 4) + (A_in - result_in(ii, 2)) *  (result_in(ii + 1, 4) - result_in(ii, 4)) / &
                                   (result_in(ii + 1, 2) - result_in(ii, 2))
                R_out = result_in(ii, 5) + (A_in - result_in(ii, 2)) * (result_in(ii + 1, 5) - result_in(ii, 5)) / &
                                   (result_in(ii + 1, 2) - result_in(ii, 2))
                exit
            endif
        enddo
    endsubroutine find_hydraulic_properties

endprogram main
