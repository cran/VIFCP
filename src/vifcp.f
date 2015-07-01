C   Main program
C      PROGRAM TEST
      subroutine vif(y,n,l,siglev,cpset,m)
C  Input: 
C  y: dataset; l:length of segment; siglev: sign. level 
      IMPLICIT NONE
      real*8 y(*),z(10000),v(10000),vtemp,vnew,c(10000)
      real*8 t,pn,rho,sigma,sigmap,sigmam
      integer i,j,k,l,n,aa,s(10000),m,temp
      integer flag,shift,length,cpset(100)
      real*8 siglev,w,dw,alpha
      integer ws
      
C   Initiate parameter values
      aa=n/l-1;      s(1)=0
      w=0.05;      dw=0.05;      flag=0
      sigmap=0.0
      sigmam=0.0 
      i=1;m=0
      do while(i<aa.and.w>0)
          alpha=w/(1+i-flag)
C   Update v(j), c(j) and vnew
          s(m+2)=i+1;
          if(m>1) then
          j=m
          else
          j=1
          end if
          do while(j<=m+1)
          k=s(j)*l+1;   vtemp=0
          do while(k<=s(j+1)*l)
          vtemp=vtemp+y(k);   k=k+1
          end do
          v(j)=vtemp;   c(j)=1.0/(s(j+1)-s(j)); j=j+1
          end do
          k=i*l+1;  vnew=0
          do while(k<=(i+1)*l)
          vnew=vnew+y(k);   k=k+1
          end do

C   Update rho and sigma
          rho=sqrt(l*(1-c(m+1)))
          if(i<2) then
          sigmap=0; k=1;
          else
          k=i*l+1;
          end if
          do while(k<=(i+1)*l)
          sigmap=sigmap+y(k)*y(k);   k=k+1
          end do
          sigmam=0; j=1
          do while(j<=m+1)
          sigmam=sigmam+c(j)*v(j)*v(j)/l; j=j+1
          end do
          sigma=sqrt((sigmap-sigmam)/((i+1)*l-m-2))
C   Compute t
          t=abs(vnew-c(m+1)*v(m+1))/rho/sigma
C   Compute pnorm(t)
          call pnorm(t,pn)
C          print *,'i=',i,', rho=', rho,', sigma=',sigma,', t=',t
C          print *,'pvalue=',pn,', alpha=',alpha

          if(2.0*pn>2.0-alpha) then
C   Test for a changepoint in (il-l/2,il+l) using the CUSUM
          shift=i*l-l/2;    length=3*l/2; k=1
          do while(k<=length)
          z(k)=y(k+shift);  k=k+1
          end do
          call changepoint(z,length,siglev,ws)
          temp=m+1;
C   If the test is significant, obtain the exact changepoint
          if(ws>0 .and. m==0) then
                  ws=ws+shift; m=m+1;  s(m+1)=i;  
                  flag=i; w=w+dw
                  cpset(m)=ws
                  temp=m;
          else if(ws>0 .and. abs(ws+shift-cpset(temp))>l)then
                  ws=ws+shift; m=m+1;  s(m+1)=i;  
                  flag=i; w=w+dw
                  cpset(m)=ws
                  temp=m;
          end if
          
          else
          ws=0
          end if
C   If there is no changepoint, the wealth will be reduced
          if(ws==0) then
          w=w-alpha/(1.0-alpha)
          end if
C   Update i until i>=aa or w<=0
          i=i+1
       end do
       return 
       end




      subroutine pnorm(t,pn)
      integer loop
      real*8 t,pn,stepsize,step
      loop=16;  stepsize=t/(2.0*loop); step=stepsize; pn=stepsize;
      do while(step<t-stepsize)
      pn=pn+4.0*exp(-step*step/2)*stepsize
      step=step+stepsize
      pn=pn+2.0*exp(-step*step/2)*stepsize
      step=step+stepsize
      end do
      pn=pn+4.0*exp(-step*step/2)*stepsize+exp(-t*t/2)*stepsize
      pn=pn/3.0/sqrt(2.0*3.1415926)+0.5
      return
      end


      subroutine changepointnew(z,n,siglev,cp)
      integer n,k,cp
      real*8 z(n),siglev,sumn,sumsquaren,sumk
      real*8 uk,ukt,umax,pv,pvta,pvtb,pvb
      k=0;  sumn=0;  sumsquaren=0
      do while(k<n)
      k=k+1;    sumn=sumn+z(k);   sumsquaren=sumsquaren+z(k)*z(k)
      end do
      k=1;  sumk=z(k);  umax=0
      do while(k<n)
      ukt=abs(n*1.0*sumk-k*1.0*sumn)
      ukb=k*(n-k)*sumsquaren-(n-k)*sumk*sumk-k*(sumn-sumk)*(sumn-sumk)
      uk=ukt/sqrt(ukb)
      if(uk>umax) then
      umax=uk;  cp=k
      end if
      k=k+1;    sumk=sumk+z(k)
      end do
      pvta=-log(-0.5*log(1-siglev))+2*log(log(n*1.0))
      pvtb=0.5*log(log(log(n*1.0)))-0.5*log(3.1415926)
      pvb=SQRT(2*log(log(n*1.0)))
      pv=(pvta+pvtb)/pvb
      if(umax<pv) then
      cp=0
      end if
      RETURN
      END


C weighted statistical testing for a signal c.p.
      subroutine changepoint (A, n, siglev,ws)

C A:        the list of data based on which to find c.p.
C n:        length of the data
C siglev:  signifiance level
      integer n
      real*8 A(n), siglev,pi
C add some useful parameters.
      real*8 Sa, bara2, bara
      real*8 ca(n-1),wa(n-1),cw(n-1)
      Integer i
      real*8 sta
      integer WS

      Sa=0.0
      bara=0.0
      bara2=0.0
      sta=0.0
      pi=3.1415926

C firstly calculate the mean of A and A^2
      do 10 i=1,n
           bara  = bara+A(i)
           bara2 = bara2+(A(i))**2
10    continue
      bara=bara/n
      bara2=bara2/n

      do 20 i=1,(n-1)
         Sa=Sa+A(i)
         ca(i)=(sa-i*bara)/(SQRT(i*(n-i)/(n*1.0)))
         wa(i)=SQRT(bara2-(Sa**2)/(i*n)-(bara*n-Sa)**2/(n*(n-i)*1.0))
         cw(i)=ABS(ca(i)/wa(i))
20    continue

      sta=cw(1)
      WS=1
      do 30 i=2,(n-1)
           IF(sta .LT. cw(i)) then 
           sta=cw(i)
           WS=i
           end if
30    continue

      pvalue=(-log(-0.5*log(1-siglev))+2*log(log(n*1.0))+
     * 0.5*log(log(log(n*1.0)))-0.5*log(pi))/
     *(SQRT(2*log(log(n*1.0))))
      IF(sta .LT. pvalue) WS=0
      RETURN
      END


C...........................................................
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
