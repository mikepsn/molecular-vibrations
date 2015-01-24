c-------------------------------------------------------------
c PROGRAM 	: ROOT FINDING
c-------------------------------------------------------------

     program main
c		call false_position(10, 1.0d0, 0.001d0)
            call new_false_position(0.33d0, 0.0d0, 0.00000001d0)
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
c-------------------------------------------------------------

      double precision function simpson(n, a, b)
      implicit none
      integer*4 n, i, factor
      real*8 a, b, answer
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
c-------------------------------------------------------------

      double precision function f(n,a,b)
      implicit none
      integer*4 n
      real*8 a, b, simpson, g
      f = simpson(n,a,b) - g(b) 
      return
      end

c-------------------------------------------------------------
c FUNCTION  : false_position
c-------------------------------------------------------------

      subroutine false_position(n,start,tolerance)  
      implicit none
      integer*4 n
      real*8 start, tolerance
      real*8 f
      real*8 h, a, b, x1, x2, x3, f1, f2, f3

      a = 0.0d0
      b = start
      h = (b-a)/dfloat(n)
      x1 = start
      x2 = x1 + h
      f1 = f(n,a, x1)
      f2 = f(n,a, x2)

      write(*,*) 'f1 = ', f1
      write(*,*) 'f2 = ', f2

      do while (f1*f2 .ge. 0.0d0) 
        x2 = x1 + h 
        f1 = f(n,a,x1) 
        f1 = f(n,a,x2) 
        write(*,*) 'f1 = ', f1
c		write(*,*) 'f2 = ', f2
      end do

      end

c-------------------------------------------------------------
c FUNCTION  : new_false_position
c-------------------------------------------------------------
      
      double precision function func(x)
      implicit none
      real*8 x
      func = x**2 - 4
      return
      end

      subroutine new_false_position(n, x1, tolx)
      implicit none
      real*8 n, x1, tolx
      real*8 func
      real*8 x2, f1, f2, x3, f3

      x2 = x1 + n
      f1 = func(x1)
      f2 = func(x2)

      do while (f1*f2 .ge. 0.0d0)
            x2 = x2 + n
            f2 = func(x2)
      end do

      x3 = x2 - f2*(x2-x1)/(f2-f1)
      f3 = func(x3)

      do while (dabs(f3) > tolx)
            if (f1*f3 .lt. 0.0d0) then
                  x2 = x3
            else
                  x1 = x3
            endif
            f1 = func(x1)
            f2 = func(x2)
            x3 = x2 - f2*(x2-x1)/(f2-f1)
            f3 = func(x3)
      end do

      write(*,*) 'x3 = ',x3

      end
      



