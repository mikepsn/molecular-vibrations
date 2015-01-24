c-------------------------------------------------------------
c PROGRAM 	: ROOT FINDING - FALSE POSITION METHOD
c 
c SUBJECT   : 640-364 COMPUTATIONAL PHYSICS
c NAME      : MICHAEL PAPASIMEON
c DATE      : 06/08/1997
c
c PURPOSE   : This program does a numerical integration 
c           : between 0 and 1 for the function ln(1+x)/x
c           : using both the Trapezoidal and Simpsons
c           : rules for quadrature.
c-------------------------------------------------------------

      program main
      implicit none
      integer*4 n, i
      real*8 result, false_position, tolx, answer, error

      answer = dsqrt(3.0d0)

      n = 10
      tolx = 1.0d0/dfloat(n)
      do 10 i = 0, 6, +1
            result = false_position(n, 1.0d0, tolx)
            error = dabs(result - answer)
            write(*,15)n, tolx, result, error
15          format(i10, f14.10, f14.10, g14.6)
            n = n*10
            tolx = 1.0d0/dfloat(n)
10    continue
      end

c-------------------------------------------------------------
c FUNCTION	: g
c PURPOSE	: returns t**2
c-------------------------------------------------------------

      double precision function g(t)
      implicit none
      real*8 t
      g = t**2
      return
      end

c-------------------------------------------------------------
c SUBROUTINE	: simpson
c PURPOSE		: calculates the integral between a and b
c               : for the function g(x)
c-------------------------------------------------------------

      double precision function simpson(n, a, b)
      implicit none
      integer*4 n, i, factor
      real*8 a, b
      real*8 h, sum, x, g

      i = 0
      factor = 4
      sum = 0.0d0
      x = 0.0d0
      h = (b-a)/dfloat(n)

      do 300 i = 1, (n-1), +1
            x = a+i*h
            if (factor .eq. 2) then
                  sum = sum + 2.0d0*g(x)
                  factor = 4
            else
                  sum = sum + 4.0d0*g(x)
                  factor = 2
            endif
300   continue
      sum = sum + g(a) + g(b)
      sum = (h*sum)/3.0d0

     simpson = sum
     return

     end

c-------------------------------------------------------------
c FUNCTION  : f
c PURPOSE   : calculates and returns the value of the function
c           : of the integral between a and b of the function
c           : g(x) determined by the simpson function minus
c           : c. When a = 0 and b = x, we have
c           : f = Integrate[t**2 dt,0,x] - x
c-------------------------------------------------------------

      double precision function f(n,a,b)
      implicit none
      integer*4 n
      real*8 a, b, simpson
      f = simpson(n,a,b) - b 
      return
      end

c-------------------------------------------------------------
c FUNCTION  : false_position
c PURPOSE   : calculates and returns the root of the function
c           : f, after position x1, to a tolerance of tolx
c           : using the false position method. The parameter
c           : n determines the accuracy to which the function
c           : f can be calculated to.
c-------------------------------------------------------------

      double precision function false_position(n, x1, tolx)
      implicit none
      integer*4 n
      real*8 x1, tolx
      real*8 f
      real*8 x2, f1, f2, x3, f3, h, a

      a = 0.0d0 
      h = 0.33

      x2 = x1 + h
      f1 = f(n,a,x1)
      f2 = f(n,a,x2)

      do while (f1*f2 .ge. 0.0d0)
            x2 = x2 + h
            f2 = f(n,a,x2)
      end do

      x3 = x2 - f2*(x2-x1)/(f2-f1)
      f3 = f(n,a,x3)

      do while (dabs(f3) .gt. tolx)
            if (f1*f3 .lt. 0.0d0) then
                  x2 = x3
            else
                  x1 = x3
            endif
            f1 = f(n,a,x1)
            f2 = f(n,a,x2)
            x3 = x2 - f2*(x2-x1)/(f2-f1)
            f3 = f(n,a,x3)
      end do

      false_position = x3
      return

      end
      

