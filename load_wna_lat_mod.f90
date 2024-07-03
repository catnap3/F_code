module load_wna_lat_mod
    implicit none
    double precision, allocatable :: thetaL_ray(:), Lat_ray(:)
    character(len=*), parameter :: FILEPATH = "../../raytracing-app/data/wave_normal_angle05.dat"
    character(len=*), parameter :: SAVEFILE = "wna_lat.dat"
    contains
    subroutine load_wna_lat
        implicit none
        double precision :: xx, yy, xd, yd, uu, vv
        integer(kind=4) :: iostat, n
		integer(kind=4) :: unit = 17
        character(len=200) :: header
        
        ! ファイル行数を取得
        open(unit=unit, file=FILEPATH, status="old", iostat=iostat)

        if (iostat /= 0) then
            print *, "Error opening file(1)."
            stop
        endif

        read(unit,*) header
        n = 0
        do
            read(unit,*, iostat=iostat)
            if (iostat /= 0) exit
            n = n + 1
        end do
        close(unit)

        ! 配列を割り当て
        allocate(thetaL_ray(n), Lat_ray(n))

        ! ファイルを再オープンしてデータを読み込み
        open(unit=unit, file=FILEPATH, status="old", iostat=iostat)
        if (iostat /= 0) then
            print *, "Error opening file(2)."
            stop
        endif

        read(unit,*) header
        n = 0
        do
            read(unit,*, iostat=iostat) xx, yy, xd, yd, uu, vv, thetaL_ray(n+1), Lat_ray(n+1)
            if (iostat /= 0) exit
            n = n + 1
        end do

        close(unit)
    end subroutine load_wna_lat

    subroutine save_data
        implicit none
        integer :: unit = 18, i
        open(unit=unit, file=SAVEFILE, status="replace")
        do i = 1, size(thetaL_ray)
            write(unit, *) thetaL_ray(i), Lat_ray(i)
        end do
        close(unit)
    end subroutine save_data

    subroutine load_data
        implicit none
        integer :: unit = 18, i, n, iostat

        open(unit=unit, file=SAVEFILE, status="old", iostat=iostat)
        if (iostat /= 0) then
            print *, "Error opening saved data file."
            stop
        endif

        n = 0
        do
            read(unit,*, iostat=iostat)
            if (iostat /= 0) exit
            n = n + 1
        end do
        close(unit)

        allocate(thetaL_ray(n), Lat_ray(n))

        open(unit=unit, file=SAVEFILE, status="old")
        do i = 1, n
            read(unit, *) thetaL_ray(i), Lat_ray(i)
        end do
        close(unit)
    end subroutine load_data
end module load_wna_lat_mod
