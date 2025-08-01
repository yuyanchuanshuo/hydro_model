program main
    implicit none
    integer :: i, number, t_all, j, dt, t,hour
    real, dimension(:, :), allocatable :: A, Q, Q0, A0, H, R
    real, dimension(:), allocatable :: n, dx, ZB
    integer :: file_unit  ! 添加文件单元号变量
    integer :: file_unit_all    ! 添加专门输出第2000时刻所有断面的文件

    number = 100 !假设有100个断面，待修改
    t_all = 100000 !计算总时段
    dt = 3 ! 计算间隔10s

    ! 分配数组内存
    allocate(A(number, t_all), A0(number, t_all), Q(number, t_all), Q0(number, t_all), R(number, t_all), H(number, t_all))!定义以断面为行，时间为列的矩阵
    allocate(n(number), dx(number), ZB(number))

    Q0(:, 1) = 100 !假设每一个断面的第一时刻流量均为100m3/s
    Q(:, 1) = 100 !假设每一个断面的第一时刻流量均为100m3/s
    A0(:, 1) = 200 !假设每一个断面的第一时刻面积均为1000m2
    A(:, 1) = 200 !假设每一个断面的第一时刻面积均为1000m2
    A(number,:) = 100 ! 假设最后一个断面的所有时刻过水面积都是200

    n = 0.035 ! 假设每一个断面的糙率均为0.035
    dx = 100 !假设每一个断面的间距为100m

    do i = 1, t_all ! 定义上边界入流，这个以后从文件中读取
        Q0(1, i) = 10 + (30 - 10) * real(i - 1) / (t_all - 1)
        Q(1, i) = 10 + (30 - 10) * real(i - 1) / (t_all - 1)
    enddo

    do i = 1, number ! 定义不同断面的初始数据，这个以后从文件中读取
        ZB(i) = 25 - (25 - 20) * real(i - 1) / (number-1)
    enddo


    ! 打开输出文件
    open(newunit=file_unit, file='output_results.txt', status='replace', action='write')
    write(file_unit, '(A)') '时间(小时)  断面编号  流量(m3/s)  面积(m2)  水深(m)'

    open(newunit=file_unit_all, file='output_all_sections_at_t2000.txt', status='replace', action='write')
    write(file_unit_all, '(A)') '断面编号  流量(m3/s)  面积(m2)  水深(m)'

    j = 88
    do t = 1, t_all - 1
        call dis_cal(number, n, dx, dt, Q0(:, t), ZB, A(:, t), Q(:, t + 1), H(:, t))  ! 使用call调用子例程
        Q0(2:number, t + 1) = Q(2:number, t + 1) ! 定义这个断面输入的Q
        Q(1,t+1) = Q0(1,t+1)
        A(:, t + 1) = area_cal(number, dt, Q(:, t+1), A0(:, t), dx)
        A0(:, t + 1) = A(:, t + 1) ! 定义这个断面输入的A

        if (mod(t,1200)==0) then
            hour = t/1200
            print *, '第', hour, '个小时第', j, '个断面的流量为：', Q(j, t), A(j, t), H(j, t)
            ! 写入文件
            write(file_unit, '(I6, I10, 3F15.6)') hour, j, Q(j, t), A(j, t), H(j, t)
        endif

        ! 在t=2000时刻输出所有断面结果
        if (t == 100000-50) then
            do i = 1, number
                write(file_unit_all, '(I10, 3F15.6)') i, Q(i, t), A(i, t), H(i, t)
            enddo
        endif
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

    subroutine dis_cal(n_size, n_in, dx_in, dt_in, Q_in, ZB_in, A_in, Q_out, h_out) !函数只能返回一个结果变量，改用子例程
        implicit none
        integer :: iii
        integer :: B = 50 ! 假设一个何宽，做试算，后期要删掉
        integer, intent(in) :: n_size, dt_in
        real, intent(in), dimension(n_size) :: dx_in, n_in
        real, dimension(n_size) :: uu, Q_in, H_in, ZB_in, A_in, R_in
        real, dimension(n_size) :: Q_out, h_out
        real :: XDX, HDX, FRI, alfa = 1 !alfa动量修正系数，一般取1

        H_in = A_in / B ! 假设一个何宽，做试算，后期要删掉
        h_out = H_in
        R_in = A_in / (2 * H_in + B) ! 假设一个何宽，做试算，后期要删掉，水力半径 = 过水面积/湿周
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

    endsubroutine dis_cal ! 子例程不需要返回值，所以去掉result


endprogram main
