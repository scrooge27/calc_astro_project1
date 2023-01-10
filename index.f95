module utils
    implicit none
    contains

    subroutine lettura(doc,matrix,n,m,position,vec)
        integer::m !m conta già le sole colonne numeriche
        integer::n !n conta il numero di righe
        real*8,allocatable::matrix(:,:)
        character(len=*)::doc
        character(len=*)::position
        character(len=*),optional,allocatable::vec(:)
        integer::i

        open(10,file=doc)
        n=0
        do
            read(10,*,end=100)
            n=n+1
        end do
        100 continue

        allocate(matrix(n,m))   
        
        if (present(vec)) then 
            allocate(vec(n))
        end if

        rewind(10)

        do  i=1,n
            select case (position)
            case("first")
                read(10,*) vec(i),matrix(i,1:m)
            case("last")
                read(10,*) matrix(i,1:m),vec(i)
            case("none")
                read(10,*) matrix(i,1:m)
            end select
        end do

        rewind(10)
        close(10)
    end subroutine lettura

    subroutine sort(v1,rows,cols,v2)
        integer::i,j,indice,rows,cols
        real*8::v1(rows,cols),min,copy1
        character(len=7),optional::v2(rows)
        character(len=7)::copy2

        
        do i=1,rows
            min=v1(i,1)
            indice=i
    
            !ricerca del minimo
            do j=i+1,rows
                if (v1(j,1)<min) then
                    min=v1(j,1)
                    indice=j
                end if
            end do
    
            v1(indice,1)=v1(i,1)      !metto il dato corrente in corrisondenza dell'indice in cui si trovava il minimo
            v1(i,1)=min              !metto nella posizione corrente il valore minimo
    
            do j=2,cols
                copy1= v1(indice,j)   !prima non avevo bisogno della copia avendo il dato salvato in min
                v1(indice,j)=v1(i,j)  !metto ogni dato in posizione corrente in corrispondenza dell'indice in cui si trovava il minimo
                v1(i,j)=copy1         !metto nella posizione corrente il dato corrispondente al valore minimo
            end do
    
            if(present(v2)) then
                copy2= v2(indice)   !prima non avevo bisogno della copia avendo il dato salvato in min
                v2(indice)=v2(i)  !metto ogni dato in posizione corrente in corrispondenza dell'indice in cui si trovava il minimo
                v2(i)=copy2     
            end if
    
        end do

    end subroutine sort

    subroutine min_fun(v,n,min,id)    !trova il minimo all'interno di un vettore
        real*8::v(n),min
        integer::n,i
        integer,optional::id
    
        min=v(1)
        if(present(id))then 
            id=1
        end if 
        do i=2,n
            if (v(i)<min) then
                min=v(i)
                if(present(id))then
                    id=i
                end if
            end if
        end do
    end subroutine min_fun

end module utils

module interpol
    implicit none
    contains
    !nelle seguenti subroutine a è la matrice dei coefficienti, c il vettore dei termini noti

    subroutine jordan(a,c,n,x)
        real*8::a(n,n),c(n),fakt,aux  
        integer::n,i,j,k
        real*8,intent(out)::x(n)
        
        do i=1,n                    !anche l'ultima variabile va normalizzata
            !normalizzazione
            aux=a(i,i)
            do j=1,n
                a(i,j)=a(i,j)/aux
            end do
            c(i)=c(i)/aux
            do j=1,n 
            !eliminazione
                if (i/=j) then
                    fakt=a(j,i)/a(i,i)
                    do k=1,n 
                        a(j,k)=a(j,k)-a(i,k)*fakt
                    end do
                    c(j)=c(j)-c(i)*fakt
                end if
            end do
        end do
        !"sostituzione"
        x=c
    end subroutine jordan

    subroutine linfit(x,y,ndati,coeff,nd)
        integer::nd, ndati
        integer::i,j,k
        real*8:: x(ndati),y(ndati)
        real*8:: a(nd,nd), c(nd), coeff(nd)
        real*8::sum
        
        do i=1,nd
           do j=1,nd
              sum=0.d0
              do k=1,ndati
                 sum=sum+x(k)**(i+j-2)
              end do
              a(i,j)=sum
            end do
        end do
      
        do i=1,nd
           sum=0.d0
           do k=1,ndati
              sum=sum+y(k)*(x(k)**(i-1))
           end do
           c(i)=sum
        end do      
        call jordan(a,c,nd,coeff)
        
    end subroutine linfit

    subroutine spline(trid,c,x,f,n)
        integer::n,i,j
        real*8::e(n-1),r(n-1),g(n-1)
        real*8::x(n),f(n)
        real*8,intent(out)::trid(n,n),c(n)
        
        c=0.d0
        trid=0.d0
        trid(1,1)=1.d0
        trid(n,n)=1.d0
        do i=2,n-1
            do j=1,n
                if(j==i)then

                    trid(i,j)=2.d0*(x(i+1)-x(i-1))

                    r(i)=trid(i,j)

                else if (j==i+1) then

                    trid(i,j)=x(i+1)-x(i)

                    g(i)=trid(i,j)
                    
                else if (j==i-1) then

                    trid(i,j)=x(i)-x(i-1)

                    e(i)=trid(i,j)
                    
                end if
            end do
            c(i)=6.d0*((f(i+1)-f(i))/(x(i+1)-x(i)))+6.d0*((f(i-1)-f(i))/(x(i)-x(i-1)))
        end do
    end subroutine spline

    real*8 function fun_spline(x,f,fd2,xval,n,i)
        integer::n,i
        real*8::xval,x(n),f(n),fd2(n)
        fun_spline= fd2(i-1)/(6.d0*(x(i)-x(i-1)))*(x(i)-xval)**3+&
                    fd2(i)/(6.d0*(x(i)-x(i-1)))*(xval-x(i-1))**3+&
                    (f(i-1)/((x(i)-x(i-1)))-fd2(i-1)*(x(i)-x(i-1))/6.d0)*(x(i)-xval)+&
                    (f(i)/((x(i)-x(i-1)))-fd2(i)*(x(i)-x(i-1))/6.d0)*(xval-x(i-1))
    end function fun_spline

end module interpol

module steps
    use utils
    use interpol
    implicit none
    contains

    subroutine step1(coeff)
        integer::n_ceph

        real*8::coeff(2)
        real*8,allocatable::mat_ceph(:,:)
        real*8,external::fun_spline
    
        character(len=17)::doc1
        character(len=7),allocatable::ceph(:)

        doc1="ceph_catalog.txt"

        call lettura(doc1,mat_ceph,n_ceph,2,"last",ceph)

        call sort(mat_ceph,n_ceph,2,ceph)

        mat_ceph(:,1)=log10(mat_ceph(:,1))
        
        call linfit(mat_ceph(:,1),mat_ceph(:,2),n_ceph,coeff,2)

        print*,coeff(1),coeff(2)
    end subroutine step1

    subroutine step2(p_vec,app_mag_vec,n_gal)
        integer::i,j,k,l,ll,ii,n_gal,n_ceph,n_p,p_id,counter,counter2,dim,dim2

        real*8,allocatable,intent(out)::p_vec(:,:),app_mag_vec(:,:)
        real*8,allocatable::mat_gal(:,:),mat_ceph(:,:),trid(:,:),diff(:),p(:),c(:),fd2(:),chisq(:),mag_arr(:),spline_vec(:)
        real*8::max_jd,min_jd,a,b,min_diff,xval(2),splval(2),min_chisq,magval,passo,jd,jd_o

        character(len=28)::spl_out
        character(len=17)::doc1,ceph
        character(len=7),allocatable::gal(:)


        doc1="gal_vel.txt"

        call lettura(doc1,mat_gal,n_gal,2,"first",gal)

        !call sort(mat_gal,n_gal,2,gal)

        allocate(p_vec(n_gal,3),app_mag_vec(n_gal,3))

        ceph(1:5)="ceph_"
        ceph(14:17)=".txt"
        
        spl_out(1:11)="spline_out_"

        do i=1,n_gal
            ceph(6:12)=gal(i)
            ceph(13:13)="_"
            !print*,"setting ",ceph
            !read(*,*)
            do j=1,3
                select case(j)
                case(1)
                    ceph(13:13)="A"
                    !print*,"ceph A"
                case(2)
                    ceph(13:13)="B"
                    !print*,"ceph B"
                case(3)
                    ceph(13:13)="C"
                    !print*,"ceph C"
                end select

                call lettura(ceph,mat_ceph,n_ceph,3,"none")
                call sort(mat_ceph,n_ceph,3)

                min_jd=mat_ceph(1,1)
                max_jd=mat_ceph(n_ceph,1)

                !carico diff con le differenze temporali
                allocate(diff(n_ceph-1))
                do k=2,n_ceph
                    diff(k-1)=mat_ceph(k,1)-mat_ceph(k-1,1)
                end do

                b=(max_jd-min_jd)*0.5d0

                !print*,"gli estremi temporali sono"
                !print*,min_jd,max_jd,b
                !read(*,*)
                
                call min_fun(diff,n_ceph-1,min_diff)
                a=3.d0*min_diff

                !carico p con il range di periodi
                n_p=int((b-a)/0.02d0)+1
                allocate (p(n_p))
                do k=1,n_p 
                    p(k)=a+(k-1)*0.02d0
                end do

                !print*,"il range di periodi da testare e il seguente"
                !print*,p(1),p(n_p),n_p
                !read(*,*)

                !traccio la spline
                allocate(trid(n_ceph,n_ceph),c(n_ceph),fd2(n_ceph),chisq(n_p))
                call spline(trid,c,mat_ceph(:,1),mat_ceph(:,2),n_ceph) !metto in relazione periodi(1) e magnitudini(2)
                call jordan(trid,c,n_ceph,fd2)

                !interpolo con la spline ogni punto traslandolo di +- P
                do k=1,n_p
                    chisq(k)=0.d0
                    counter=0
                    do l=1,n_ceph
                        xval(1)=mat_ceph(l,1)+p(k)
                        xval(2)=mat_ceph(l,1)-p(k)

                        do ll=1,2

                            do while ((xval(ll)<=max_jd).and.(xval(ll)>=min_jd))
                                
                                do ii=2,n_ceph
                                    if((xval(ll)<=mat_ceph(ii,1)).and.(xval(ll)>=mat_ceph(ii-1,1))) then 
                                        exit
                                    end if
                                end do
                                splval(ll)=fun_spline(mat_ceph(:,1),mat_ceph(:,2),fd2,xval(ll),n_ceph,ii)
                            
                                !calcolo del chi quadro
                            
                                chisq(k)=chisq(k)+((splval(ll)-mat_ceph(l,2))/mat_ceph(l,3))**2

                                select case (ll)
                                case(1)
                                    xval(ll)=xval(ll)+p(k)
                                case(2)
                                    xval(ll)=xval(ll)-p(k)
                                end select

                                counter=counter+1

                            end do
                        end do
                    end do
                    chisq(k)=chisq(k)/(counter-1)         !ho tolto il vincolo del valore stimato con la spline       
                end do

                call min_fun(chisq,n_p,min_chisq,p_id)
                !print*,"il chi quadro minimo e",min_chisq
                !print*,"il periodo minimo e",p(p_id)
                !read(*,*)

                p_vec(i,j)=p(p_id)

                !calcolo delle magnitudini 
                counter=0
                passo=p(p_id)*0.01d0
                dim=int(b*2.d0/passo)+1
                allocate(spline_vec(dim))

                dim2=int(b*2.d0/p(p_id))
                allocate(mag_arr(dim2))
                

                jd=min_jd
                jd_o=min_jd
                magval=0.d0
                counter2=2
                l=1
                k=1
                spl_out(12:28)=ceph
                OPEN(10,file=spl_out)
                
                do while(jd<=max_jd)

                    do ii=counter2,n_ceph
                        if((jd<=mat_ceph(ii,1)).and.(jd>=mat_ceph(ii-1,1))) then 
                            exit
                        end if
                    end do

                    if(ii>counter2)then
                        counter2=ii
                    end if
                    
                    spline_vec(k)=fun_spline(mat_ceph(:,1),mat_ceph(:,2),fd2,jd,n_ceph,ii)
                    WRITE(10,*)jd,spline_vec(k)
                    magval=magval+spline_vec(k)
                    k=k+1
                    jd=jd+passo
                    counter=counter+1
                    if (jd-jd_o>=p(p_id)) then
                        jd_o=jd
                        mag_arr(l)=magval/counter  !carico la media delle magnitudini nel range coperto
                        magval=0.d0
                        counter=0
                        l=l+1
                    end if
                end do

                close(10)

                magval=0.d0

                do k=1,dim2
                    magval=magval+mag_arr(k)
                end do
                
                magval=magval/dim2

                app_mag_vec(i,j)=magval

                !print*,p_vec(i,j),app_mag_vec(i,j)
                !read(*,*)
                !last operation
                deallocate(p,diff,fd2,c,trid,chisq,mat_ceph,mag_arr,spline_vec)
                
            end do
        end do
    end subroutine step2

end module steps

program main
    use steps
    implicit none
    integer::n_gal,i,j
    real*8::coeff(2)
    real*8,allocatable::p(:,:),app_mag(:,:)

    call step1(coeff)
    read(*,*)
    call step2(p,app_mag,n_gal)
    do i=1,n_gal
        print*,"galaxy number",i
        do j=1,3
            print*,p(i,j),app_mag(i,j)
        end do
        read(*,*)
    end do
    
end program main


