     Program tsit10
      real*8 w0, pi, xend, h, w(50),c(50),d(50),q(50), p(50)
      real*8 dVdq(5), dTdp(50)
      integer i,j,k,isteps,n
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  INTEGRATION of 2-BODY
c  
c  q    : displacements
c  p    : velocities
c  T    : kinetic energy (here T=.5*p'.p)
c  V    : potential (here V=-1/||q||)
c  H    : Hamiltonian, H=T(p)+V(q)
c  
c  at each step we
c  integrate from time x (where p_0,q_0) to time x+h (for p_34,q_34)
c  using
c  q_{i+1}=q_i+c_i*h*dT/dp
c  p_{i+1}=p_i-d_i*h*dV/dq
c
c  here we use eccentricity 0.5, so initialy q_0=[.5 0] and p_0=[0 sqrt(3)]
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c the coefficients
      w(1)=0.02690013604768968151437422144441685467297755661d0
      w(2)=0.939801567135683337900037741418097674632505563d0
      w(3)=-0.00803583920385358749646880826318191806393661063d0
      w(4)=-0.866485197373761372803661767454208401679117010d0
      w(5)=0.1023112911193598731078563285067131541328142449d0
      w(6)=-0.1970772151393080101376018465105491958660525085d0
      w(7)=0.617877713318069357335731125307019691019646679d0
      w(8)=0.1907272896000121001605903836891198270441436012d0
      w(9)=0.2072605028852482559382954630002620777969060377d0
      w(10)=-0.395006197760920667393122535979679328161187572d0
      w(11)=-0.582423447311644594710573905438945739845940956d0
      w(12)=0.742673314357319863476853599632017373530365297d0
      w(13)=0.1643375495204672910151440244080443210570501579d0
      w(14)=-0.615116639060545182658778437156157368647972997d0
      w(15)=0.2017504140367640350582861633379013481712172488d0
      w(16)=0.45238717224346720617588658607423353932336395045d0
      w0=1.0d0-2.0d0*(w(1)+w(2)+w(3)+w(4)+w(5)+w(6)+w(7)+w(8)+w(9)+
     1w(10)+w(11)+w(12)+w(13)+w(14)+w(15)+w(16))
      c(1)=0.5d0*w(16)
      c(34)=0.5d0*w(16)
      c(2)=0.5d0*(w(16)+w(15))
      c(33)=0.5d0*(w(16)+w(15))
      c(3)=0.5d0*(w(15)+w(14))
      c(32)=0.5d0*(w(15)+w(14))
      c(4)=0.5d0*(w(14)+w(13))
      c(31)=0.5d0*(w(14)+w(13))
      c(5)=0.5d0*(w(13)+w(12))
      c(30)=0.5d0*(w(13)+w(12))
      c(6)=0.5d0*(w(12)+w(11))
      c(29)=0.5d0*(w(12)+w(11))
      c(7)=0.5d0*(w(11)+w(10))
      c(28)=0.5d0*(w(11)+w(10))
      c(8)=0.5d0*(w(10)+w(9))
      c(27)=0.5d0*(w(10)+w(9))
      c(9)=0.5d0*(w(9)+w(8))
      c(26)=0.5d0*(w(9)+w(8))
      c(10)=0.5d0*(w(8)+w(7))
      c(25)=0.5d0*(w(8)+w(7))
      c(11)=0.5d0*(w(7)+w(6))
      c(24)=0.5d0*(w(7)+w(6))
      c(12)=0.5d0*(w(6)+w(5))
      c(23)=0.5d0*(w(6)+w(5))
      c(13)=0.5d0*(w(5)+w(4))
      c(22)=0.5d0*(w(5)+w(4))
      c(14)=0.5d0*(w(4)+w(3))
      c(21)=0.5d0*(w(4)+w(3))
      c(15)=0.5d0*(w(3)+w(2))
      c(20)=0.5d0*(w(3)+w(2))
      c(16)=0.5d0*(w(2)+w(1))
      c(19)=0.5d0*(w(2)+w(1))
      c(17)=0.5d0*(w(1)+w0)
      c(18)=0.5d0*(w(1)+w0)
      d(34)=0.0d0
      d(33)=w(16)
      d(1)=w(16)
      d(32)=w(15)
      d(2)=w(15)
      d(31)=w(14)
      d(3)=w(14)
      d(30)=w(13)
      d(4)=w(13)
      d(29)=w(12)
      d(5)=w(12)
      d(28)=w(11)
      d(6)=w(11)
      d(27)=w(10)
      d(7)=w(10)
      d(26)=w(9)
      d(8)=w(9)
      d(25)=w(8)
      d(9)=w(8)
      d(24)=w(7)
      d(10)=w(7)
      d(23)=w(6)
      d(11)=w(6)
      d(22)=w(5)
      d(12)=w(5)
      d(21)=w(4)
      d(13)=w(4)
      d(20)=w(3)
      d(14)=w(3)
      d(19)=w(2)
      d(15)=w(2)
      d(18)=w(1)
      d(16)=w(1)
      d(17)=w0
c
c the driver
c    initial values
       pi=4.d0*datan(1.d0)
       xend=10*pi
       isteps=155
       n=2
       q(1)=.5d0
       q(n)=0.d0
       p(1)=0.d0
       p(n)=dsqrt(3.d0)
       h=xend/isteps
c
c    the main loop
        do 1821 i=1,isteps
           do 1822 j=1,34
              call kinet(q,p,n,dTdp)
              do 1823 k=1,n
                 q(k)=q(k)+c(j)*h*dTdp(k)
 1823         continue
              call poten(q,p,n,dVdq)
              do 1824 k=1,n
                 p(k)=p(k)-d(j)*h*dVdq(k)
 1824         continue
 1822      continue
 1821 continue
      write(*,*) q(1),q(2),p(1),p(2)
      stop
      end
c
      subroutine kinet(q,p,n,dTdp)
      integer n
      real*8 q(50), p(50), dTdp(50)
      do 17 i=1,n
  17  dTdp(i)=p(i)
      return
      end

      subroutine poten(q,p,n,dVdq)
      integer n
      real*8 q(50), p(50), dVdq(50)
      do 18 i=1,n
  18  dVdq(i)=q(i)/dsqrt(q(1)**2+q(2)**2)**3
      return
      end
