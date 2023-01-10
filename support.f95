!salvo in un array una serie di punti interpolati per fare il grafico
                dim=int(b*200.d0)
                allocate(spline_vec(dim))
                counter2=2
                jd=min_jd

                spl_out(12:28)=ceph
                OPEN(10,file=spl_out)

                do l=1,dim
                    do ii=counter2,n_ceph
                        if((jd<=mat_ceph(ii,1)).and.(jd>=mat_ceph(ii-1,1))) then 
                            exit
                        end if
                    end do

                    if(ii>counter2)then
                        counter2=ii 
                    end if

                    spline_vec(l)=fun_spline(mat_ceph(:,1),mat_ceph(:,2),fd2,jd,n_ceph,ii)
                    WRITE(10,*)jd,spline_vec(l)

                    jd=jd+0.01d0
                end do
                close(10)