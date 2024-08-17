program main 
implicit none
real :: a1(1:3), a2(1:3), a3(1:3), x(1:3), shell
real ::  dr
real, allocatable :: g(:,:), rij(:,:,:), ri(:,:)
real norm,volume
integer :: i, j, ct, k1, k2, k3
integer :: elenum(1:2)
character(len=2)  elename(1:2)
real :: pi
parameter(pi = 3.1415926535)       

open(unit = 226, file = 'POSCAR_a-SiO2')
read(226,*)
read(226,*)
read(226,*) a1
read(226,*) a2
read(226,*) a3
read(226,*)elename
read(226,*)elenum
allocate(g(int(norm(a1+a2+a3)/0.016)+2,2))
allocate(ri(sum(elenum),3))
allocate(rij(sum(elenum),sum(elenum),3))
read(226,*)
read(226,*)

do  i = 1, int(norm(a1+a2+a3)/0.016)+2
    dr = -0.008 + 0.016*i
    g(i,1) = dr
end do

do  i = 1, sum(elenum)
    read(226,*) x
    do j = 1,3
        if (x(j)>= 0.5) then
            x(j) = x(j) - 1.0
        end if
    end do
    ri(i,:) = x
end do

do i = 1,sum(elenum)
    do j = i+1, sum(elenum)
        rij(i,j,:) =    (ri(j,1)-ri(i,1))*a1 &
                       +(ri(j,2)-ri(i,2))*a2 &
                       +(ri(j,3)-ri(i,3))*a3
    end do
end do

open(unit=2271, file='RDF_SiO2_OO')
do ct = 1, int(norm(a1+a2+a3)/0.016)+2
    g(ct,2) = 0.0
    do i = 1,elenum(1)-1
        do j = i+1,elenum(1)
            do k1 = -2,2 
                do k2 = -2,2
                    do k3 = -2,2
            x = rij(i,j,:) + k1*a1 + k2*a2 + k3*a3
            if (abs(g(ct,1)-norm(x)) <= 0.008) then
                g(ct,2) = g(ct,2) + 1.0
            end if
                    end do
                end do
            end do
        end do
    end do
!   print *, g(ct,:)
    g(ct,2) = g(ct,2)/elenum(1)
    shell   = (4*pi/3)*(0.008)*( 3*g(ct,1)**2 + 0.008**2/4)
    g(ct,2) = (g(ct,2)/shell)/(elenum(1)/volume(a1,a2,a3))
    write(2271,*) g(ct,:)
    print *, real(ct)/(int(norm(a1+a2+a3)/0.016)+2)*100,'%'
end do
print *,'O-O is done!'
close(2271)

open(unit=2272, file='RDF_SiO2_SiSi')
do ct = 1, int(norm(a1+a2+a3)/0.016)+2
    g(ct,2) = 0.0
    do i = elenum(1)+1, elenum(1)+elenum(2)-1
        do j = i+1, elenum(1)+elenum(2)
            do k1 = -2,2 
                do k2 = -2,2
                    do k3 = -2,2
            x = rij(i,j,:) + k1*a1 + k2*a2 + k3*a3
            if (abs(g(ct,1)-norm(x)) <= 0.008) then
                g(ct,2) = g(ct,2) + 1.0
            end if
                    end do
                end do
            end do
        end do
    end do
!   print *, g(ct,:)
    g(ct,2) = g(ct,2)/elenum(2)
    shell   = (4*pi/3)*(0.008)*( 3*g(ct,1)**2 + 0.008**2/4)
    g(ct,2) = (g(ct,2)/shell)/(elenum(2)/volume(a1,a2,a3))
    write(2272,*) g(ct,:)
    print *, real(ct)/(int(norm(a1+a2+a3)/0.016)+2)*100,'%'
end do
print *,'Si-Si is done!'
close(2272)

open(unit=2272, file='RDF_SiO2_SiO')
do ct = 1, int(norm(a1+a2+a3)/0.016)+2
    g(ct,2) = 0.0
    do i = 1, elenum(1)
        do j = elenum(1)+1, elenum(1)+elenum(2)
            do k1 = -2,2 
                do k2 = -2,2
                    do k3 = -2,2
            x = rij(i,j,:) + k1*a1 + k2*a2 + k3*a3
            if (abs(g(ct,1)-norm(x)) <= 0.008) then
                g(ct,2) = g(ct,2) + 1.0
            end if
                    end do
                end do
            end do
        end do
    end do
!   print *, g(ct,:)
    g(ct,2) = g(ct,2)/elenum(1)
    shell   = (4*pi/3)*(0.008)*( 3*g(ct,1)**2 + 0.008**2/4)
    g(ct,2) = (g(ct,2)/shell)/(elenum(1)/volume(a1,a2,a3))
    write(2272,*) g(ct,:)
    print *, real(ct)/(int(norm(a1+a2+a3)/0.016)+2)*100,'%'
end do
print *,'Si-O is done!'
close(2272)
print *, 'OK!'
end

function volume(v1,v2,v3)
real, intent(in) :: v1(1:3), v2(1:3), v3(1:3)
real :: volume
volume = v1(1)*v2(2)*v3(3) + v1(2)*v2(3)*v3(1) + v1(3)*v2(1)*v3(2) &
       - v1(1)*v2(3)*v3(2) - v1(2)*v2(1)*v3(3) - v1(3)*v2(2)*v3(1)
return
end function volume

function norm(v)
real, intent(in) :: v(1:3)
real :: norm
norm = (v(1)**2 + v(2)**2 + v(3)**2)**0.5
return
end function norm
