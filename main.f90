! Truong Dang
! 01-11-2016
program Laplace_2D
    implicit none
    interface
       function functionf2D(x,y,k)
          real(8)::functionf2D,y,x
          integer(4)::k
       end function functionf2D

       function exact_solution2D(x,y,k)
          real(8)::exact_solution2D,y,x
          integer(4)::k
       end function exact_solution2D
    end interface

    integer(4)::N,k1,j,i,t
    real(8)::aa,bb
    real(8)::h,k
    real(8),dimension(4)::M
    real(8),dimension(:),allocatable::x,xcp,y,ycp,F,u
    real(8),dimension(:,:),allocatable::B,A,C,D,B_temp,B_temp_inv
    real(8),dimension(:,:),allocatable::u_dis,u_exact

    real(8),dimension(4)::norm_l2,norm_l1

    real(8)::ai,bi,cj,dj,sij
    integer(4)::i1,i2

    N=2
    k1=3

    do j=1,4
       aa=0.0d0
       bb=1.0d0
       N=N*2
       allocate(x(1:N))
       allocate(xcp(1:N+1))
       allocate(y(1:N))
       allocate(ycp(1:N+1))
       allocate(B(1:(N-1)*(N-1),1:(N-1)*(N-1)))
       allocate(B_temp(1:(N-1)*(N-1),1:(N-1)*(N-1)))
       allocate(B_temp_inv(1:(N-1)*(N-1),1:(N-1)*(N-1)))
       allocate(A(1:N-1,1:N-1))
       allocate(C(1:N-1,1:N-1))
       allocate(D(1:N-1,1:N-1))
       allocate(F(1:(N-1)*(N-1)))
       allocate(u(1:(N-1)*(N-1)))
       allocate(u_dis(1:N+1,1:N+1))
       allocate(u_exact(1:N+1,1:N+1))

       M(j)=N

    !%%Make mesh points x(i+1/2)
    h=1.0/(N-1)
    x(1)=aa
    do i=2,N-1
        x(i)=(i-1)*h + x(1)
    enddo
    x(N)=bb

    !%%Make control points xcp(i)
    xcp=0.0d0
    xcp(1)=x(1)
    do i=2,N+1
        if (i==N+1) then
            xcp(i)=x(i-1)
        else
            xcp(i)=(x(i)+x(i-1))/2.0d0;
        endif
    enddo

    !%%Make mesh points y(i+1/2) i=0,1,...,N
    k=1.0d0/(N-1)
    y(1)=aa
    do i=2,N-1
        y(i)=(i-1)*k + y(1)
    enddo
    y(N)=bb

    !%%Make control points ycp(i) i=0,1,...,N+1
    ycp=0.0d0
    ycp(1)=y(1)
    do i=2,N+1
        if (i==N+1) then
            ycp(i)=y(i-1)
        else
            ycp(i)=(y(i)+y(i-1))/2.0d0
        endif
    enddo
    !%%

    B=0.0d0



    !%%
    ai=-1.0d0/h**2
    bi=-1.0d0/h**2
    cj=-1.0d0/k**2
    dj=-1.0d0/k**2
    sij=(ai+bi+cj+dj)
    !%%
    A=0.0d0
    do i=1,N-1
        if(i==1) then
            A(i,i)=sij
            A(i,i+1)=-bi
            elseif(i==N-1)then
                A(i,i)=sij
                A(i,i-1)=-ai
            else
                A(i,i)=sij
                A(i,i-1)=-ai
                A(i,i+1)=-bi
            endif

    enddo
    !%%

        C=0.0d0
        do i=1,N-1
            C(i,i)=-cj
        enddo

    !%%

        D=0.0d0
        do i=1,N-1
            D(i,i)=-dj
        enddo

    !%%

    do i=1,N-1
        if(i==1)then
            B((i-1)*(N-1)+1:i*(N-1),(i-1)*(N-1)+1:i*(N-1))=-A
            B((i-1)*(N-1)+1:i*(N-1),i*(N-1)+1:(i+1)*(N-1))=-D
        elseif(i==N-1)then
                B((i-1)*(N-1)+1:i*(N-1),(i-1)*(N-1)+1:i*(N-1))=-A
                B((i-1)*(N-1)+1:i*(N-1),(i-2)*(N-1)+1:(i-1)*(N-1))=-C
            else
                B((i-1)*(N-1)+1:i*(N-1),(i-1)*(N-1)+1:i*(N-1))=-A
                B((i-1)*(N-1)+1:i*(N-1),i*(N-1)+1:(i+1)*(N-1))=-D
                B((i-1)*(N-1)+1:i*(N-1),(i-2)*(N-1)+1:(i-1)*(N-1))=-C

        endif
    enddo
    !%A=(1/h^2)*A;
    F=0.0d0
    do i1=1,N-1
        do i2=1,N-1
            F((i1-1)*(N-1)+i2)=functionf2D(xcp(i2),ycp(i1),k1)
        enddo
    enddo
    !u=B\F;
    B_temp = B
    call gauss_2(B,F,u,(N-1)*(N-1))
    call matinv(B_temp,(N-1)*(N-1),B_temp_inv)
    !print*,B_temp_inv

    !do i1=1,(N-1)*(N-1)
    !        u(i1) = (1.0D0/B(i1,i1))*(F(i1) &
    !        -dot_product(B(i1,1:(N-1)*(N-1)),u)+B(i1,i1) &
    !        *u(i1))
    !enddo


    u_dis=0.0d0
    do i1=1,N-1
        do i2=1,N-1
            u_dis(i1+1,i2+1)=u((i1-1)*(N-1)+i2)
        enddo
    enddo

    open(unit=1,file='u_dis')
    do i1=1,N+1
        do i2=1,N+1
            write(1,*),xcp(i1),ycp(i2),u_dis(i1,i2)
        enddo
        write(1,*),''
    enddo
    close(unit=1)



    u_exact=0.0d0
    do i1=1,N+1
        do i2=1,N+1
            u_exact(i1,i2)=exact_solution2D(xcp(i2),ycp(i1),k1)
        enddo
    enddo

    open(unit=2,file='u_exact')
    do i1=1,N+1
        do i2=1,N+1
            write(2,*),xcp(i1),ycp(i2),u_exact(i1,i2)
        enddo
        write(2,*),''
    enddo
    close(unit=2)
    !figure
    !surf(xcp,ycp,u_dis)
    !title('u discrete', 'FontWeight','bold')
    !if(k1==2)
    !    zlim([0 0.08])
    !else
    !    zlim([0 0.025])
    !end

    !figure
    !surf(xcp,ycp,u_exact)
    !title('u exact', 'FontWeight','bold')
    !%% Odd indices figures: u discrete
    !%% Even indices figures: u exact

    !%%
    norm_l2(j)=0.0d0
    norm_l1(j)=0.0d0

    do i=1,N+1
        do t=1,N+1
        norm_l2(j)=norm_l2(j)+(u_dis(i,t)-u_exact(i,t))**2*h*k
        norm_l1(j)=norm_l1(j)+abs(u_dis(i,t)-u_exact(i,t))*h*k
        enddo
    enddo
    norm_l2(j)=sqrt(norm_l2(j))
    print*,'Mesh size:',N,'x',N
    print*,'L2 norm = ',norm_l2(j)
    print*,'L1 norm = ',norm_l1(j)

!    %%
!    %%
!    norm_h1(j)=0;
!    %norm_h1(j)=norm_h1(j)+(((u_dis(1)-u_exact(1))))^2;
!     for i=1:N
!         for t=1:N
!             if (i==1)&&(t==1)
!                 norm_h1(j)=norm_h1(j)+4*(((u_dis(i,t)-u_exact(i,t))))^2;
!             else
!                 norm_h1(j)=norm_h1(j)+(abs(u_dis(i+1,t)-u_exact(i+1,t))-abs(u_dis(i,t)-u_exact(i,t)))^2+(abs(u_dis(i,t+1)-u_exact(i,t+1))-abs(u_dis(i,t)-u_exact(i,t)))^2;
!             end
!         end
!     end
!     %norm_h1(j)=norm_h1(j)+(((u_dis(N)-u_exact(N))-(u_dis(N-1)-u_exact(N-1))))^2;
!     norm_h1(j)=(norm_h1(j))^(1/2);
!
!    %%


!figure
!%plot( log(M), -log(norm_l2), 'red', log(M), 2*log(M),'black', log(M), -log(norm_h1), 'cyan', log(M), (1/2)*log(M), 'green')
!%ylim([0 20])
!%legend('|| ||_{L^2}','2*Log(N)','|| ||_{H^1_0}','(1/2)*Log(N)')
!%legend('Location','northwest')
!plot( log(M), -log(norm_l2))
!hold on
!plot( log(M), 2*log(M),'r')
!xlabel("Log(N)")
!ylabel("-Log(err)")
    deallocate(x)
    deallocate(xcp)
    deallocate(y)
    deallocate(ycp)
    deallocate(B)
    deallocate(B_temp)
    deallocate(B_temp_inv)
    deallocate(A)
    deallocate(C)
    deallocate(D)
    deallocate(F)
    deallocate(u)
    deallocate(u_dis)
    deallocate(u_exact)
enddo

open(unit=3,file='L2 norm')
do j=1,4
    write(3,*),log(M(j)),log(norm_l2(j))
enddo
close(unit=3)

open(unit=4,file='L1 norm')
do j=1,4
    write(4,*),log(M(j)),log(norm_l1(j))
enddo
close(unit=4)


end program

function functionf2D(x,y,k)
   implicit none
   real(8)::functionf2D,y,x
   integer(4)::k

if (k==1) then
    functionf2D =2
endif

if (k==2) then
    functionf2D =2*y*(1-y)+2*x*(1-x)
endif

if (k==3) then
    functionf2D =2*x**2*(2*y - 2)*(x - 1) + 2*x**2*y*(x - 1) + 4*x*y*(y - 1)**2 + 2*y*(x - 1)*(y - 1)**2
endif

!if (k==4)
!    f=2*y^2*(2*x - 2)*(y - 1) + 4*x*y*(x - 1)^2 + 2*x*y^2*(y - 1) + 2*x*(x - 1)^2*(y - 1);
!end

if (k==5) then
    functionf2D =0.0d0
endif

end function functionf2D

function exact_solution2D(x,y,k)
    implicit none
    real(8)::exact_solution2D,y,x
    integer(4)::k
if (k==1) then
    exact_solution2D =x*(1-x)
endif
!% Truong added 1-11-2016
!% Neu truyen k=1, ham chi phu thuoc x nen nghiem xap xi ra sai

if(k==2)then
    exact_solution2D=x*(1-x)*y*(1-y)
endif

if (k==3) then
    exact_solution2D =x**2*(1-x)*y*(1-y)**2
endif

!if (k==4) then
!    u_ex=x*(1-x)^2*y^2*(1-y)
!endif

if (k==5) then
    exact_solution2D = sin(x)*exp(y)
endif
end function exact_solution2D


