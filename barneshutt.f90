!****************************************************************************************
!
!  Created by Michael Larson
!  Feburary 1st 2020
!
!
!****************************************************************************************


!****************************************************************************************
! 
!****************************************************************************************
module barnes_hutt
  use, intrinsic :: iso_c_binding
  implicit none

  type :: Vec3
     real(8),dimension(1:3) :: v
  end type Vec3
  
  type :: SimData
     integer :: count
     type(Vec3),pointer,dimension(:) :: pos, vel, acc
     real(8),pointer,dimension(:) :: mass
  end type SimData

  type :: BHNode
     integer :: count
     integer,pointer,dimension(:) :: indices
     type(Vec3) :: min,max
     type(Vec3) :: pos,vel,acc
     type(BHNode),pointer :: parent
     type(BHNode),pointer,dimension(:) :: nodes
  end type BHNode
  
  interface
     function drand48() bind(C,name="drand48")
       use, intrinsic :: iso_c_binding
       real(c_double) :: drand48
     end function drand48
  end interface

contains
  subroutine setVec3(v,x,y,z)
    type(Vec3) :: v
    real(8) :: x, y, z

    v%v(1) = x
    v%v(2) = y
    v%v(3) = z
  end subroutine setVec3

  function makeVec3(x, y, z)
    real(8) :: x, y, z
    type(Vec3) :: makeVec3

    makeVec3%v(1) = x
    makeVec3%v(2) = y
    makeVec3%v(3) = z    
  end function makeVec3
  
  function vec3Add(a, b)
    type(Vec3) :: a, b
    type(Vec3) :: vec3Add

    vec3Add%v = a%v + b%v
  end function vec3Add
  
  function vec3Sub(a, b)
    type(Vec3) :: a, b
    type(Vec3) :: vec3Sub

    vec3Sub%v = a%v - b%v
  end function vec3Sub

  function vec3Mul(a, b)
    type(Vec3) :: a
    real(8) :: b
    type(Vec3) :: vec3Mul

    vec3Mul%v = a%v * b
  end function vec3Mul
 
  function drand48Vec3()
    type(Vec3) :: drand48Vec3

    drand48Vec3%v(1) = drand48()
    drand48Vec3%v(2) = drand48()
    drand48Vec3%v(3) = drand48()
  end function drand48Vec3
    
  subroutine initData(data, count)
    type(SimData) :: data
    integer(c_int) :: count

    integer :: i 
    
    data%count = count

    allocate(data%pos(1:count))
    allocate(data%vel(1:count))
    allocate(data%acc(1:count))
    allocate(data%mass(1:count))

    do i=1,count
       data%pos(i) = drand48Vec3()
       data%vel(i) = drand48Vec3()
       data%acc(i) = drand48Vec3()
       data%mass(i) = drand48()
    end do

  end subroutine initData

  subroutine sort_bodies(bh_node, data)
    type(BHNode) :: bh_node
    type(SimData) :: data

    integer :: i, j, index, count, bin, node_count
    type(Vec3) :: min, max, mid, delta, temp_vec3
    real(8) :: one = 1.0D0
    
    count = bh_node%count

    ! allocate sub-nodes of this node along with enough for
    ! all indices
    allocate(bh_node%nodes(1:8))

    do i=1,8
       bh_node%nodes(i)%count = 0
       allocate(bh_node%nodes(i)%indices(1:count))
    end do
    
    ! figure out mid point for bh_node sorting
    min = bh_node%min
    max = bh_node%max
    delta = vec3Sub(max, min)
    delta = vec3Mul(delta, 0.5D0)
    mid = vec3Add(min, delta)

    ! init sub node min and maxes
    !
    !max -----------
    !    |    |    |
    !    |    |    |
    !mid -----------
    !    |    |    |
    !    |    |    |
    !min -----------
    !   min mid  max
    
    bh_node%nodes(1)%min = makeVec3(min%v(1), min%v(2), min%v(3))
    bh_node%nodes(2)%min = makeVec3(mid%v(1), min%v(2), min%v(3))

    bh_node%nodes(3)%min = makeVec3(min%v(1), mid%v(2), min%v(3))
    bh_node%nodes(4)%min = makeVec3(mid%v(1), mid%v(2), min%v(3))

    bh_node%nodes(5)%min = makeVec3(min%v(1), min%v(2), mid%v(3))
    bh_node%nodes(6)%min = makeVec3(mid%v(1), min%v(2), mid%v(3))
    
    bh_node%nodes(7)%min = makeVec3(min%v(1), mid%v(2), mid%v(3))
    bh_node%nodes(8)%min = makeVec3(mid%v(1), mid%v(2), mid%v(3))

    ! create the node max by adding mid to min for that node
    do i=1,8
       bh_node%nodes(i)%max = vec3Add(bh_node%nodes(i)%min, delta)
    end do

    ! sort each body to a node
    do i=1,count
       ! indices store the index into data
       index = bh_node%indices(i)
       temp_vec3 = data%pos(index)

       bin = 1
       
       if (temp_vec3%v(1) > mid%v(1)) bin = bin + 1
       
       if (temp_vec3%v(2) > mid%v(2)) bin = bin + 2
       
       if (temp_vec3%v(3) > mid%v(3)) bin = bin + 4

       node_count = bh_node%nodes(bin)%count + 1
       
       bh_node%nodes(bin)%indices(node_count) = index

       bh_node%nodes(bin)%count = node_count
    end do
  end subroutine sort_bodies

end module barnes_hutt


program main
  use, intrinsic :: iso_c_binding
  use barnes_hutt
  implicit none

  type(SimData),pointer :: data
  type(BHNode),pointer :: bh_node
  type(BHNode),pointer,dimension(:,:) :: queue, temp_queue
  integer :: i, j, k, count, q_size, q_index, q_b1, q_b2, q_entries, new_q_entries
  real(8) :: zero = 0.0D0, one = 1.0D0
  
  allocate(data)
  allocate(bh_node)
  
  count = 1000000
  call initData(data, count)

  bh_node%count = count
  allocate(bh_node%indices(1:count))

  do i=1,count
     bh_node%indices(i) = i
  end do

  bh_node%min = makeVec3(zero, zero, zero)
  bh_node%max = makeVec3(one, one, one)

  ! this is a double buffered queue
  ! entires are added to one buffer then
  ! entires are submitted from the other buffer
  ! then the buffer is swapped...
  q_size = 256
  allocate(queue(2,q_size))

  q_index = 1
  q_b1 = 1
  q_b2 = 2
  queue(q_b1, q_index) = bh_node
  do
     !$OMP PARALLEL DO
     submit_enrites:do i=1,q_index
        call sort_bodies(queue(q_b1, i), data)
     end do submit_enrites
     !$OMP END PARALLEL DO

     q_entries = q_index
     q_index = 0
     
     !print *,"q entries",q_entries,"q_size",q_size
 
     each_entry:do i=1,q_entries
        each_node:do j=1,8
           ! add sub-nodes to queue for next round
           add_entry:if (queue(q_b1, i)%nodes(j)%count > 1) then
              q_index = q_index + 1
              queue(q_b2, q_index) = queue(q_b1, i)%nodes(j)
              
              !print *,"index",i,"add",q_index,"count",queue(q_b1, i)%nodes(j)%count
              !print "(f6.3 f6.3 f6.3)",queue(q_b1, i)%nodes(j)%min
              
              ! reallocate queue if queue grows too large
              full_queue:if (q_index >= q_size) then
                 !print *, "full queue, size", q_size

                 ! copy current queue to temp queue
                 allocate(temp_queue(1:2, q_size))

                 do k=1,q_entries
                    temp_queue(q_b1,k) = queue(q_b1,k)
                 end do

                 do k=1,q_index
                    temp_queue(q_b2,k) = queue(q_b2,k)
                 end do

                 deallocate(queue)

                 ! resize queue and copy back from temp queue
                 q_size = q_size * 2
                 allocate(queue(1:2, q_size))

                 do k=1,q_entries
                    queue(q_b1,k) = temp_queue(q_b1,k)
                 end do

                 do k=1,q_index
                    queue(q_b2,k) = temp_queue(q_b2,k)
                 end do

                 deallocate(temp_queue)
              end if full_queue
           end if add_entry
        end do each_node
     end do each_entry

     !print *, "q_index", q_index

     ! no nodes with more than 1 left
     if (q_index == 0) exit

     !print *, "swapping queue buffers", q_b1, q_b2

     ! swap q buffers
     if (q_b1 == 1) then
        q_b1 = 2
        q_b2 = 1
     else
        q_b1 = 1
        q_b2 = 2
     end if
  end do
     
end program main
