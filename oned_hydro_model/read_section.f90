program read_section
    implicit none
    integer :: file_unit, file_unit_out  ! 添加文件单元号变量
    integer :: n_sections = 0, n_points = 0, ios, i, j, section_num
    integer :: total_lines = 0 ! 文件总行数
    integer, dimension(:), allocatable :: id_line_numbers  ! 新增：存储ID行行号
    character(len=256) :: line, start_distance_str

    type :: section_data
        character(len=256) :: id
        real :: start_distance
        real, dimension(:), allocatable :: distances
        real, dimension(:), allocatable :: elevations

    endtype section_data

    type(section_data), dimension(:), allocatable :: sections

    ! 在原有type定义后添加新的数据类型
    type :: hydraulic_data
        character(len=256) :: section_id
        real, dimension(:,:), allocatable :: results  ! 存储水位、面积、宽度、湿周、水力半径
    end type hydraulic_data

     ! 在变量声明部分添加
    type(hydraulic_data), dimension(:), allocatable :: section_hydraulics


    ! 打开文件并读取数据
    open(newunit=file_unit, file='data_input/cross_section_data.csv', status='old', action='read')

    ! 第一遍读取确定断面数量
    do
        read(file_unit, '(A)', iostat=ios) line
        if(ios /= 0) exit
        total_lines = total_lines + 1
    enddo
    rewind(file_unit)

    ! 分配ID行行号数组
    allocate(id_line_numbers(total_lines))  ! 最大可能数量
    n_sections = 0

    do i = 1, total_lines
        read(file_unit, '(A)') line
        ! 判断是否是断面ID行(行首是字符串)
        if(verify(line(1:1), '0123456789') > 0) then  ! 简化判断方法
            n_sections = n_sections + 1
            id_line_numbers(n_sections) = i
        endif
    enddo
    ! 压缩数组，移除未使用的0元素
    if(n_sections < total_lines) then
        id_line_numbers = pack(id_line_numbers, id_line_numbers /= 0)
    endif

    rewind(file_unit) ! 回到文件起点

    ! 分配断面数据内存
    allocate(sections(n_sections))
    section_num = 0
    j = 0
    ! 第二遍读取：逐行读取，判断是ID行，n就+1
    do i = 1, total_lines
        if(i == id_line_numbers(section_num + 1)) then
            if(i /= 1) then
                ! 调整数组大小为实际读取的点数
                sections(section_num)%distances = sections(section_num)%distances(1:n_points)
                sections(section_num)%elevations = sections(section_num)%elevations(1:n_points)
            endif

            section_num = section_num + 1
            ! 解析断面ID和起点距
            read(file_unit, '(A)', iostat=ios) line  ! 读取整行到line变量
            read(line, *) sections(section_num)%id, start_distance_str! 按空格分隔读取多个变量
            sections(section_num)%start_distance = parse_distance(start_distance_str)
            n_points = 0
            j = j + 1
        else
            ! 确保在使用前分配了数组内存
            if(.not. allocated(sections(section_num)%distances)) then
                allocate(sections(section_num)%distances(1000))
                allocate(sections(section_num)%elevations(1000))
            endif
            n_points = n_points + 1
            read(file_unit, '(A)', iostat=ios) line  ! 读取整行到line变量
            read(line, *) sections(section_num)%distances(n_points), sections(section_num)%elevations(n_points)
            j = j + 1

            ! 最后一行读取后调整section的大小
            if(i == total_lines) then
                ! 调整数组大小为实际读取的点数
                sections(section_num)%distances = sections(section_num)%distances(1:n_points)
                sections(section_num)%elevations = sections(section_num)%elevations(1:n_points)
            endif
        endif

    enddo
    close(file_unit)

    ! 添加断面倒序判断逻辑，确保最后一个断面最低
    if (minval(sections(1)%elevations) < minval(sections(n_sections)%elevations)) then
        sections = sections(n_sections:1:-1)  ! 倒序断面数组
    endif


    ! 新增：输出断面详细数据到文件
    open(newunit=file_unit, file='data_output/section_details.txt', status='replace', action='write')
    write(file_unit, '(A)') '断面序号 断面ID 起点距(m) 点号 距离(m) 高程(m)'
    
    do i = 1, n_sections
        do j = 1, size(sections(i)%distances)
            write(file_unit, '(I10, A10, F15.3, I10, 2F15.3)') i, trim(sections(i)%id), &
                sections(i)%start_distance, &
                j, &
                sections(i)%distances(j), &
                sections(i)%elevations(j)
        enddo
    enddo
    close(file_unit)


    ! 在断面倒序判断后添加
    allocate(section_hydraulics(n_sections))

    ! 修改计算和存储部分
    do i = 1, n_sections
        section_hydraulics(i)%section_id = sections(i)%id
        call calculate_hydraulic_properties(sections(i), section_hydraulics(i)%results)
    enddo

    ! 新增：输出到文件
    open(newunit=file_unit_out, file='data_output/section_hydraulic_properties.txt', status='replace', action='write')
    write(file_unit_out, '(A)') '断面ID 水位(m) 面积(m2) 水面宽度(m) 湿周(m) 水力半径(m)'
    
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

contains
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

    ! 新增函数：计算断面水力特性
    subroutine calculate_hydraulic_properties(section, results_matrix)
        type(section_data), intent(in) :: section
        real, dimension(:, :), allocatable, intent(out) :: results_matrix
        real :: min_elev, max_elev, current_h
        integer :: iii, jjj, n_levels, n_points_in
        real :: left_dist, right_dist, left_elev, right_elev
        real :: segment_width, segment_area, segment_perimeter
        real, dimension(:), allocatable :: areas, widths, wetted_perimeters
        
        n_points_in = size(section%distances)
        min_elev = minval(section%elevations)
        max_elev = maxval(section%elevations)
        
        ! 设置水位范围：从最低点到最高点+10m，步长0.1m
        n_levels = int((max_elev + 10.0 - min_elev) / 0.1) + 1
        allocate(results_matrix(n_levels, 5))  ! 4列：水位、面积、宽度、湿周、水力半径
        allocate(areas(n_levels))
        allocate(widths(n_levels))
        allocate(wetted_perimeters(n_levels))
        
        ! 初始化水位数组
        do iii = 1, n_levels
            results_matrix(iii, 1) = min_elev + (iii-1)*0.1
        enddo
        
        ! 计算每个水位对应的水力特性
        do iii = 1, n_levels
            current_h = results_matrix(iii, 1)
            areas(iii) = 0.0
            widths(iii) = 0.0
            wetted_perimeters(iii) = 0.0
            
            ! 遍历所有断面点
            do jjj = 2, n_points_in
                left_dist = section%distances(jjj-1)
                right_dist = section%distances(jjj)
                left_elev = section%elevations(jjj-1)
                right_elev = section%elevations(jjj)
                
                ! 计算当前水位下的有效宽度
                if (current_h > min(left_elev, right_elev)) then
                    segment_width = right_dist - left_dist   

                    ! 计算梯形面积
                    if (current_h >= max(left_elev, right_elev)) then
                        segment_area = segment_width * (current_h - (left_elev + right_elev)/2.0)
                        segment_perimeter = segment_width + &
                            sqrt((right_dist - left_dist)**2 + (right_elev - left_elev)**2)
                    else
                        ! 部分淹没情况(三角形)
                        if (left_elev < right_elev) then
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
            if (wetted_perimeters(iii) > 0.0) then
                results_matrix(iii, 5) = areas(iii) / wetted_perimeters(iii)
            else
                results_matrix(iii, 5) = 0.0
            endif
        enddo
        
        deallocate(areas, widths, wetted_perimeters)
    end subroutine calculate_hydraulic_properties

endprogram
