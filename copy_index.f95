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

    subroutine max_fun(v,n,max,id)    !trova il massimo all'interno di un vettore
        real*8::v(n),max
        integer::n,i
        integer,optional::id
    
        max=v(1)
        if(present(id))then 
            id=1
        end if 
        do i=2,n
            if (v(i)>max) then
                max=v(i)
                if(present(id))then
                    id=i
                end if
            end if
        end do

    end subroutine max_fun

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

    subroutine closest(x,y,n,id,xvals)          
        real*8::x(n),y(n),xvals(2),yvals(2),seek
        real*8,parameter::passo=1.d0
        integer::n,id,i,min_id(2)
        
        seek=y(id)+passo
        
        call max_fun(y,n,yvals(1)) !in questo modo alla prima iterazione il secondo termine dell'and risulta sempre vero
        yvals(2)=yvals(1)

        !cerco a sx
        do i=1,id-1
            if(abs(y(i)-seek)<=abs(yvals(1)-seek))then
                yvals(1)=y(i)
                min_id(1)=i
            end if 
        end do
        
        !cerco a dx
        do i=id+1,n
            if(abs(y(i)-seek)<=abs(yvals(2)-seek))then
                yvals(2)=y(i)
                min_id(2)=i
            end if 
        end do

        xvals(1)=x(min_id(1))

        xvals(2)=x(min_id(2))
    
        !print*,min_id(1),min_id(2)
    end subroutine closest

    subroutine w_mean(v,ve,n,m,me)
        integer::n,i
        real*8::v(n),ve(n),w(n),sv,sw
        real*8,intent(out)::m,me
        w=(1/ve)**2
        sv=0.d0
        sw=0.d0
        do i=1,n
            sv=sv+v(i)*w(i)
            sw=sw+w(i)
        end do
        m=sv/sw
        me=1/sqrt(sw)
    end subroutine w_mean

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

module integration
    implicit none
    contains
    subroutine gen(a)
        real*8::a,b,risultato,sum,x_int
        integer::n_punti,n_int
        real*8,external::f,f1

        sum=0.d0

        print*,"dammi il numero di punti"
        read(*,*) n_punti
        print*,"dammi il numero di intervallini"        !spezzo su più intervalli come Newton-Cotes
        read(*,*) n_int
        print*,"dammi il valore dove spezzare l'integrale"
        read(*,*) x_int

        !faccio l'integrale con l'infinito negativo
        a=-1.d0/abs(x_int)
        b=0.d0

        call quadratura(f1,a,b,n_punti,n_int,risultato)

        print*,risultato

        sum=sum+risultato

        !faccio l'integrale centrale
        a=-abs(x_int)
        b=abs(x_int)

        call quadratura(f,a,b,n_punti,n_int,risultato)

        print*,risultato

        sum=sum+risultato

        !faccio l'integrale con l'infinito positivo
        a=0.
        b=1./abs(x_int)
        

        call quadratura(f1,a,b,n_punti,n_int,risultato)
        
        print*,risultato

        sum=sum+risultato

        print*,sum

    end subroutine gen

    subroutine quadratura(f,a,b,n,n_int,res)
        real*8,external::f
        real*8,allocatable::c(:),xd(:),x(:)
        real*8::a,b,res,sum,delta,x1,x2
        integer::n,n_int,i,j
    
        allocate(c(0:n-1), xd(0:n-1), x(0:n-1))     !li rendo zero-based soltanto per coerenza rispetto alle slide
    
        select case(n)
        case(1)
            xd(0)=0.d0
            c(0)=2.d0
    
        case(2)
             xd(0)=-sqrt(1.d0/3.d0)
            xd(1)=-xd(0)
            c(0)=1.d0
            c(1)=1.d0
    
        case(3)
            xd(0)=-sqrt(3.d0/5.d0)
            xd(1)=0.
            xd(2)=-xd(0)
            c(0)=5.d0/9.d0
            c(1)=8.d0/9.d0
            c(2)=c(0)
    
        case(4)
            xd(0)=-sqrt(3.d0/7.d0-2.d0/7.d0*sqrt(6.d0/5.d0))
            xd(1)=-xd(0)
            xd(2)=-sqrt(3.d0/7.d0+2.d0/7.d0*sqrt(6.d0/5.d0))
            xd(3)=-xd(2)
            c(0)=(18.d0+sqrt(30.d0))/36.d0
            c(1)=c(0)
            c(2)=(18.d0-sqrt(30.d0))/36.d0
            c(3)=c(2)
    
            
            
        case default
            print*,"case not implemented yet"
            stop
        end select
    
        delta=(b-a)/n_int
        res=0.d0
    
        do j=1,n_int
            x1=a+(j-1)*delta
            x2=a+j*delta 
    
            do i=0,n-1
                x(i)=((x2+x1)+(x2-x1)*xd(i))/2.d0
            end do
    
            sum=0.d0
    
            do i=0,n-1
                sum=sum+c(i)*f(x(i))
            end do
    
            sum=sum*(x2-x1)/2.d0  !moltiplico per dx
            
            res=res+sum
    
        end do
       
    end subroutine quadratura

end module integration

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

    end subroutine step1

    subroutine step2(ceph,coeff,p,app_mag,abs_mag,dist_mod,dist)
        integer::i,j,k,ii,n_ceph,n_p,p_id,counter,counter2,counter3,dim,dim2

        real*8,intent(in)::coeff(2)
        real*8,intent(out)::p(2),app_mag(2),abs_mag(2),dist_mod(2),dist(2)
        real*8,allocatable::mat_ceph(:,:),trid(:,:),diff(:),p_vec(:),c(:),fd2(:),chisq(:),mag_arr(:),mag_err(:),spline_vec(:)
        real*8::max_jd,min_jd,a,b,min_diff,xval(2),splval(2),min_chisq,cl_p(2),magval,passo,jd

        character(len=17),intent(in)::ceph
        character(len=23)::chisq_file
        character(len=28)::spl_out

        !inizializzo i file di output per i grafici

        chisq_file(1:6)="chisq_"
        chisq_file(7:23)=ceph

        spl_out(1:11)="spline_out_"
        spl_out(12:28)=ceph

        call lettura(ceph,mat_ceph,n_ceph,3,"none")
        call sort(mat_ceph,n_ceph,3)

        min_jd=mat_ceph(1,1)
        max_jd=mat_ceph(n_ceph,1)

        !carico diff con le differenze temporali
        allocate(diff(n_ceph-1))
        do i=2,n_ceph
            diff(i-1)=mat_ceph(i,1)-mat_ceph(i-1,1)
        end do

        b=(max_jd-min_jd)*0.5d0

        !print*,"gli estremi temporali sono"
        !print*,min_jd,max_jd,b
        !read(*,*)
                
        call min_fun(diff,n_ceph-1,min_diff)
        a=3.d0*min_diff

        !carico p_vec con il range di periodi
        n_p=int((b-a)/0.02d0)+1
        allocate (p_vec(n_p))
        do i=1,n_p 
            p_vec(i)=a+(i-1)*0.02d0
        end do

        !print*,"il range di periodi da testare e il seguente"
        !print*,p_vec(1),p_vec(n_p),n_p
        !read(*,*)

        !traccio la spline
        allocate(trid(n_ceph,n_ceph),c(n_ceph),fd2(n_ceph),chisq(n_p))
        call spline(trid,c,mat_ceph(:,1),mat_ceph(:,2),n_ceph) !metto in relazione periodi(1) e magnitudini(2)
        call jordan(trid,c,n_ceph,fd2)

        open(10,file=chisq_file)
            !interpolo con la spline ogni punto traslandolo di +- P
        do i=1,n_p
            chisq(i)=0.d0
            counter=0
            do j=1,n_ceph
                xval(1)=mat_ceph(j,1)+p_vec(i)
                xval(2)=mat_ceph(j,1)-p_vec(i)

                do k=1,2

                    do while ((xval(k)<=max_jd).and.(xval(k)>=min_jd))
                        
                        do ii=2,n_ceph
                            if((xval(k)<=mat_ceph(ii,1)).and.(xval(k)>=mat_ceph(ii-1,1))) then 
                                exit
                            end if
                        end do
                        splval(k)=fun_spline(mat_ceph(:,1),mat_ceph(:,2),fd2,xval(k),n_ceph,ii)
                    
                        !calcolo del chi quadro
                    
                        chisq(i)=chisq(i)+((splval(k)-mat_ceph(j,2))/mat_ceph(j,3))**2

                        select case (k)
                        case(1)
                            xval(k)=xval(k)+p_vec(i)
                        case(2)
                            xval(k)=xval(k)-p_vec(i)
                        end select

                        counter=counter+1

                    end do
                end do
            end do
            chisq(i)=chisq(i)/(counter-1)         !ho tolto il vincolo del valore stimato con la spline 
            write(10,*)p_vec(i),chisq(i)      
        end do

        close(10)

        call min_fun(chisq,n_p,min_chisq,p_id)
        !print*,"il chi quadro minimo e",min_chisq
        !print*,"il periodo minimo e",p_vec(p_id),"con id",p_id
        !read(*,*)
        call closest(p_vec,chisq,n_p,p_id,cl_p)
        !print*,"il periodo a sinistra è",cl_p(1)
        !print*,"il periodo a destra è",cl_p(2)
        !read(*,*)
        cl_p(1)=p_vec(p_id)-cl_p(1)
        cl_p(2)=cl_p(2)-p_vec(p_id)

        p(1)=p_vec(p_id)                !periodo
        p(2)=min(cl_p(1),cl_p(2))       !errore sul periodo


        !calcolo della magnitudine apparente
        counter=0
        passo=p(1)*0.01d0
        dim=int(b*2.d0/passo)+1
        
        allocate(spline_vec(dim))

        dim2=int((dim-1)/100)+1

        allocate(mag_arr(dim2))
        allocate(mag_err(dim2))

        mag_err=0.d0

        jd=min_jd
        magval=0.d0
        counter2=2
        counter3=0
        i=1
        j=1
        
        open(10,file=spl_out)
                
        do while(jd<=max_jd)

            do ii=counter2,n_ceph
                if((jd<=mat_ceph(ii,1)).and.(jd>=mat_ceph(ii-1,1))) then 
                    exit
                end if
            end do

            if(ii>counter2)then
                counter2=ii
            end if
            
            spline_vec(j)=fun_spline(mat_ceph(:,1),mat_ceph(:,2),fd2,jd,n_ceph,ii)
            write(10,*)jd,spline_vec(j)
            magval=magval+spline_vec(j)

            if ((jd-min_jd)>=i*p(1) .or. (jd+passo)>max_jd) then
                mag_arr(i)=magval/counter  !carico la media delle magnitudini nel range coperto
                !calcolo errore
                do k=j,counter3,-1
                    mag_err(i)=mag_err(i)+(spline_vec(k)-mag_arr(i))**2
                end do
                mag_err(i)=sqrt(mag_err(i)/(counter-1))
                magval=0.d0
                counter=0
                counter3=j+1
                i=i+1
            end if
            j=j+1
            counter=counter+1
            jd=jd+passo
        end do
        
        close(10)

        call w_mean(mag_arr,mag_err,dim2,app_mag(1),app_mag(2))
        
        !print*,p,app_mag
        !read(*,*)


        abs_mag(1)=coeff(2)*log10(p(1))+coeff(1)            !mag assoluta
        abs_mag(2)=abs(coeff(2)/(p(1)*log(10.d0)))*p(2)     !errore sulla magnitudine
        

        dist_mod(1)=app_mag(1)-abs_mag(1)                   !modulo di distanza
        dist_mod(2)=sqrt(app_mag(2)**2+abs_mag(2)**2)       !errore sul modulo di distanza


        dist(1)=10**(0.2d0*dist_mod(1)-5)                   !distanza
        dist(2)=dist(1)*log(10.d0)*0.2d0*dist_mod(2)        !errore sulla distanza

        close(10)
        !last operation
        deallocate(p_vec,diff,fd2,c,trid,chisq,mat_ceph,mag_arr,mag_err,spline_vec)
    end subroutine step2

    subroutine step5(h0)
        !da implementare
        real*8,parameter::t0=13.88d9,t0_err=0.01*t0
        real*8::h0(2)
        
    end subroutine step5
    
end module steps

program main
    use steps
    implicit none
    integer::n_gal,i,j
    real*8::coeff(2),h0(2), &                                        !variabili generali
            p(2),app_mag(2),abs_mag(2),dist_mod(3,2),dist(3,2)      !variabili delle cefeidi
    real*8,allocatable::gal_vel(:,:), &                             !variabili delle galassie
                        gal_dist(:,:), &
                        gal_h(:,:) 
    character(len=11)::doc1="gal_vel.txt"
    character(len=7),allocatable::gal(:)
    character(len=17)::ceph

    call lettura(doc1,gal_vel,n_gal,2,"first",gal)
    allocate(gal_dist(n_gal,2),gal_h(n_gal,2))
    call step1(coeff)

    print*,coeff
    !nota: c1 è c(2) e c2 è c(1)
    
    ceph(1:5)="ceph_"
    ceph(14:17)=".txt"

    open(20,file="tab_ceph.txt")
    open(15,file="tab_gal.txt")

    do i=1,n_gal
        ceph(6:12)=gal(i)
        ceph(13:13)="_"
        print*,"setting ",ceph
        read(*,*)
        do j=1,3
            select case(j)
            case(1)
                ceph(13:13)="A"
                print*,"ceph A"
            case(2)
                ceph(13:13)="B"
                print*,"ceph B"
            case(3)
                ceph(13:13)="C"
                print*,"ceph C"
            end select
            call step2(ceph,coeff,p,app_mag,abs_mag,dist_mod(j,:),dist(j,:))  
            write(20,*) p,app_mag,abs_mag,dist_mod(j,:),dist(j,:)
            print*,"period:",p(1),"error:",p(2)
            print*,"apparent magnitude:",app_mag(1),"error:",app_mag(2)     !c'è ancora un problema in questa stima!!!
            print*,"absolute magnitude:",abs_mag(1),"error:",abs_mag(2)
            print*,"distance module:",dist_mod(j,1),"error:",dist_mod(j,2)
            print*,"distance:",dist(j,1),"error:",dist(j,2)
        end do
        print*,""
        call w_mean(dist(:,1),dist(:,2),3,gal_dist(i,1),gal_dist(i,2))      !calcolo la distanza della galassia ed il suo errore
        gal_h(i,1)=gal_vel(i,1)/gal_dist(i,1)                               !calcolo la costante di Hubble per la galassia                                
        gal_h(i,2)=sqrt((gal_vel(i,2)/gal_dist(i,1))**2+ &                  !calcolo l'errore sulla costante di Hubble
                    (gal_vel(i,1)*gal_dist(i,2)/(gal_dist(i,1))**2)**2)
        write(15,*) gal_dist(i,:),gal_vel(i,:),gal_h(i,:)
        print*,"galaxy distance:",gal_dist(i,1),"error:",gal_dist(i,2)
        print*,"galaxy velocity:",gal_vel(i,1),"error:",gal_vel(i,2)
        print*,"galaxy H0:",gal_h(i,1),"error:",gal_h(i,2)
    end do  

    close(15)
    close(20)

    call w_mean(gal_h(:,1),gal_h(:,2),n_gal,h0(1),h0(2))
    print*,"mean Hubble constant:",h0(1),"error:",h0(2)

    call step5(h0)
end program main

    real*8 function f(z)
    implicit none
        real*8::z,e
        real*8,parameter:: om=0.5d0,ol=0.5d0
        e=om*(1+z)**3+(1-om-ol)*(1+z)**2+ol
        f=1/((1+z)*sqrt(e))
    end function f

    real*8 function f1(t)
    implicit none
        real*8::t
        real*8,external::f
        f1=f(1.d0/t)/t**2
    end function f1