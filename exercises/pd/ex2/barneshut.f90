module barneshut

    use particle
    use geometry

    implicit none

    real(8), parameter:: theta = 0.75

    type range 
        real(8), dimension(3) :: min, max
    end type range

    type CPtr
        type(cell), pointer :: ptr
    end type CPtr

    type cell
        type(range):: range
        type(particle3d):: part
        integer :: pos
        integer :: type !! 0 = no particle; 1 = particle; 2 = conglomerado
        real(8):: mass
        type(point3d) :: c_o_m
        type(CPtr), dimension(2,2,2):: subcell
    end type cell

    

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate_Ranges !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Calcula los rangos de las part´ıculas en la
!! matriz r en las 3 dimensiones y lo pone en la
!! variable apuntada por goal
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Calculate_Ranges(goal, part)
        type(cell), pointer :: goal
        type(particle3d), dimension(:), intent(in) :: part
        real(8), dimension(3) :: mins, maxs, medios
        real(8):: span

        mins(1) = minval(part%p%x, dim=1)
        mins(2) = minval(part%p%y, dim=1)
        mins(3) = minval(part%p%z, dim=1)

        maxs(1) = maxval(part%p%x, dim=1)
        maxs(2) = maxval(part%p%y, dim=1)
        maxs(3) = maxval(part%p%z, dim=1)
        ! Al calcular span le sumo un 10% para que las
        ! particulas no caigan justo en el borde
        span = maxval(maxs-mins) * 1.1
        medios = (maxs + mins) / 2.0
        goal%range%min = medios - span/2.0
        goal%range%max = medios + span/2.0
    end subroutine Calculate_Ranges

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Find_Cell !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Encuentra la celda donde colocaremos la particula.
!! Si la celda que estamos considerando no tiene
!! particula o tiene una particula, es esta celda donde
!! colocaremos la particula.
!! Si la celda que estamos considerando es un "conglomerado",
!! buscamos con la funci´on BELONGS a que subcelda de las 8
!! posibles pertenece y con esta subcelda llamamos de nuevo
!! a Find_Cell
!!
!! NOTA: Cuando se crea una celda "conglomerado" se crean las
!! 8 subceldas, por lo que podemos asumir que siempre existen
!! las 8. Las celdas vac´ıas se borran al final del todo, cuando
!! todo el ´arbol ha sido ya creado.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    recursive subroutine Find_Cell(root,goal,part)
        type(particle3d) :: part
        type(cell), pointer :: root, goal, temp
        integer :: i, j, k

        select case (root%type)
            case(2)
                out: do i = 1,2
                    do j = 1,2
                        do k = 1,2
                            if (Belongs(part, root%subcell(i,j,k)%ptr)) then
                                call Find_Cell(root%subcell(i,j,k)%ptr,temp,part)
                                goal => temp
                                exit out
                            end if
                        end do
                    end do
                end do out
            case default
                goal => root
        end select

    end subroutine Find_Cell
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Place_Cell !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Se ejecuta tras Find_Cell, en la celda que
!! esa funci´on nos devuelve, por lo que siempre
!! es una celda de tipo 0 (sin particula) o de tipo 1
!! (con una particula). En el caso de que es una celda
!! de tipo 1 habra que subdividir la celda y poner en
!! su lugar las dos particulas (la que originalmente
!! estaba, y la nueva).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    recursive subroutine Place_Cell(goal,part,n)
        type(cell), pointer :: goal,temp
        type(particle3d):: part
        integer:: n
    

        select case (goal%type)
            case (0)
                goal%type = 1
                goal%part = part
                goal%pos = n
            case (1)
                call Crear_Subcells(goal)
                call Find_Cell(goal,temp,part)
                call Place_Cell(temp,part,n)
            case default
                print*,"SHOULD NOT BE HERE. ERROR!"
        end select
    end subroutine Place_Cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Crear_Subcells !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Esta funcion se llama desde Place_Cell y
    !! solo se llama cuando ya hay una particula
    !! en la celda, con lo que la tenemos que
    !! subdividir. Lo que hace es crear 8 subceldas
    !! que "cuelgan" de goal y la particula que
    !! estaba en goal la pone en la subcelda que
    !! corresponda de la 8 nuevas creadas.
    !!
    !! Para crear las subceldas utilizar las funciones
    !! CALCULAR_RANGE, BELONGS y NULLIFY_POINTERS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Crear_Subcells(goal)
        type(cell), pointer :: goal
        type(particle3d):: part
        integer :: i,j,k
        integer, dimension(3) :: octant

        part = goal%part
        goal%type=2
        
        do i = 1,2
            do j = 1,2
                do k = 1,2
                    octant = (/i,j,k/)
                    allocate(goal%subcell(i,j,k)%ptr)
                    goal%subcell(i,j,k)%ptr%range%min = Calcular_Range(0,goal,octant)
                    goal%subcell(i,j,k)%ptr%range%max = Calcular_Range(1,goal,octant)
                    if (Belongs(part,goal%subcell(i,j,k)%ptr)) then
                        goal%subcell(i,j,k)%ptr%part = part
                        goal%subcell(i,j,k)%ptr%type = 1
                        goal%subcell(i,j,k)%ptr%pos = goal%pos
                    else
                        goal%subcell(i,j,k)%ptr%type = 0
                    end if
                    call Nullify_Pointers(goal%subcell(i,j,k)%ptr)
                end do
            end do
        end do
    end subroutine Crear_Subcells

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Nullify_Pointers !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Simplemente me NULLIFYca los punteros de
!! las 8 subceldas de la celda "goal"
!!
!! Se utiliza en el bucle principal y por
!! CREAR_SUBCELLS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Nullify_Pointers(goal)
        type(cell), pointer :: goal
        integer :: i,j,k

        do i = 1,2
            do j = 1,2
                do k = 1,2
                    NULLIFY(goal%subcell(i,j,k)%ptr)
                end do
            end do
        end do
    end subroutine Nullify_Pointers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Belongs !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Devuelve TRUE si la particula "part" est´a
!! dentro del rango de la celda "goal"
!!
!! Utilizada por FIND_CELL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function Belongs(part,goal)
        type(particle3d) :: part
        type(cell), pointer :: goal
        logical :: Belongs

        if (part%p%x >= goal%range%min(1) .and. &
            part%p%x < goal%range%max(1) .and. &
            part%p%y >= goal%range%min(2) .and. &
            part%p%y < goal%range%max(2) .and. &
            part%p%z >= goal%range%min(3) .and. &
            part%p%z < goal%range%max(3)) then
            Belongs = .true.
        else
            Belongs = .false.
        end if
    end function Belongs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calcular_Range !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Dado un octante "otctant" (1,1,1, 1,1,2 ... 2,2,2),
!! calcula sus rangos en base a los rangos de
!! "goal". Si "what" = 0 calcula los minimos. Si what=1
!! calcula los maximos.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function Calcular_Range(what,goal,octant)
        integer:: what
        type(cell), pointer:: goal
        integer, dimension(3) :: octant
        real(8), dimension(3):: Calcular_Range, valor_medio

        valor_medio = (goal%range%min + goal%range%max) / 2.0

        select case (what)
            case(0)
                where (octant == 1)
                    Calcular_Range = goal%range%min
                elsewhere
                    Calcular_Range = valor_medio
                end where
            case(1)
                where (octant == 1)
                    Calcular_Range = valor_medio
                elsewhere
                    Calcular_Range = goal%range%max
                end where
        end select
    end function Calcular_Range


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Borrar_empty_leaves !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Se llama una vez completado el arbol para
!! borrar (DEALLOCATE) las celdas vac´ıas (i.e.
!! sin part´ıcula).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    recursive subroutine Borrar_empty_leaves(goal)
        type(cell), pointer:: goal
        integer:: i,j,k

        if (associated(goal%subcell(1,1,1)%ptr)) then
            do i = 1,2
                do j = 1,2
                    do k = 1,2
                        call Borrar_empty_leaves(goal%subcell(i,j,k)%ptr)
                        if (goal%subcell(i,j,k)%ptr%type == 0) then
                            deallocate(goal%subcell(i,j,k)%ptr)
                        end if
                    end do
                end do
            end do
        end if
    end subroutine Borrar_empty_leaves

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Borrar_tree !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Borra el arbol completo, excepto la "head".
!!
!! El arbol se ha de regenerar continuamente,
!! por lo que tenemos que borrar el antiguo
!! para evitar "memory leaks".
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    recursive subroutine Borrar_tree(goal)
        type(cell), pointer :: goal
        integer :: i, j ,k

        do i = 1,2
            do j = 1,2
                do k = 1,2
                    if (associated(goal%subcell(i,j,k)%ptr)) then
                        call Borrar_tree(goal%subcell(i,j,k)%ptr)
                        deallocate(goal%subcell(i,j,k)%ptr)
                    end if
                end do
            end do
        end do
    end subroutine Borrar_tree

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate_masses !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Nos calcula para todas las celdas que cuelgan
!! de "goal" su masa y su center-of-mass.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    recursive subroutine Calculate_masses(goal, part)
        type(cell), pointer:: goal
        integer:: i,j,k
        real(8):: mass
        type(particle3d), dimension(:), intent(in) :: part

        goal%mass = 0
        goal%c_o_m = point3d(0,0,0)

        select case (goal%type)
            case (1)
                goal%mass = part(goal%pos)%m
                goal%c_o_m = part(goal%pos)%p
            case (2)
                do i = 1,2
                    do j = 1,2
                        do k = 1,2
                            if (associated(goal%subcell(i,j,k)%ptr)) then
                                call Calculate_masses(goal%subcell(i,j,k)%ptr, part)
                                mass = goal%mass
                                goal%mass = goal%mass + goal%subcell(i,j,k)%ptr%mass
                                goal%c_o_m%x = (mass * goal%c_o_m%x + goal%subcell(i,j,k)%ptr%mass &
                                * goal%subcell(i,j,k)%ptr%c_o_m%x) / goal%mass
                                goal%c_o_m%y = (mass * goal%c_o_m%y + goal%subcell(i,j,k)%ptr%mass &
                                * goal%subcell(i,j,k)%ptr%c_o_m%y) / goal%mass
                                goal%c_o_m%z = (mass * goal%c_o_m%z + goal%subcell(i,j,k)%ptr%mass &
                                * goal%subcell(i,j,k)%ptr%c_o_m%z) / goal%mass
                            end if
                        end do
                    end do
                end do
        end select
    end subroutine Calculate_masses

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate_forces !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Calcula las fuerzas de todas las particulas contra "head".
!! Se sirve de la funcion Calculate_forces_aux que es la
!! que en realidad hace los calculos para cada particula
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Calculate_forces(head, part)
        type(cell), pointer:: head
        integer:: i,n
        type(particle3d),  dimension(:):: part

        n = size(part)
        do i = 1,n
            call Calculate_forces_aux(i,head, part)
        end do
    end subroutine Calculate_forces


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Calculate_forces_aux !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Dada una particula "goal" calcula las fuerzas
    !! sobre ella de la celda "tree". Si "tree" es una
    !! celda que contiene una sola particula el caso
    !! es sencillo, pues se tratan de dos particulas.
    !!
    !! Si "tree" es una celda conglomerado, hay que ver primero
    !! si l/D < theta. Es decir si el lado de la celda (l)
    !! dividido entre la distancia de la particula goal
    !! al center_of_mass de la celda tree (D) es menor que theta.
    !! En caso de que asi sea, tratamos a la celda como un
    !! sola particula. En caso de que no se menor que theta,
    !! entonces tenemos que considerar todas las subceldas
    !! de tree y para cada una de ellas llamar recursivamente
    !! a Calculate_forces_aux
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    recursive subroutine Calculate_forces_aux(goal,tree, part)
        type(cell), pointer:: tree
        integer:: i,j,k, goal
        real(8):: l,D, r2, r3
        type(particle3d),   dimension(:):: part
        type(vector3d):: rji

        select case (tree%type) 
            case (1) 
                if (goal .NE. tree%pos) then
                    rji = vector_between(part(goal)%p, tree%c_o_m)
                    r2 = rji%x**2 + rji%y**2 + rji%z**2
                    r3 = r2 * sqrt(r2) + 10**(-6)
                    part(goal)%a = vsum(part(goal)%a, mulvr(part(tree%pos)%m, divvr(rji, r3)))
                end if  
            case (2)
                l = tree%range%max(1) - tree%range%min(1)
                rji = vector_between(part(goal)%p, tree%c_o_m)
                r2 = rji%x**2 + rji%y**2 + rji%z**2
                D = sqrt(r2)
                if (l/D < theta) then
                    r3 = r2 * D + 10**(-6)
                    part(goal)%a = vsum(part(goal)%a, mulvr(tree%mass, divvr(rji, r3)))
                else
                    do i = 1,2
                        do j = 1,2
                            do k = 1,2
                                if (associated(tree%subcell(i,j,k)%ptr)) then
                                    call Calculate_forces_aux(goal,tree%subcell(i,j,k)%ptr,part)
                                end if
                            end do
                        end do
                    end do
                end if
        end select
                       

    end subroutine Calculate_forces_aux
end module barneshut
