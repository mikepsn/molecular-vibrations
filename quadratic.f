c-------------------------------------------------------------
c PROGRAM	      : main
c-------------------------------------------------------------

      program main
      implicit none
      integer i, n, step
      real*8 gamma, en, answer, error
      real*8 result, tol, false_position

      write(*,*)'Enter n, gamma, answer'
      read(*,*) n, gamma, answer

      en = -0.1d0
      step = 10
      tol = 1.0d0/step

      write(*,*) '--------------------------------------------'
      write(*,*) 'Molecular Vibrations : Quadratic Potential'
      write(*,*) ' '
      write(*,*) 'step = step size for simpson integration'
      write(*,*) 'tol = tolerance for false position method'
      write(*,*) 'n = ',n
      write(*,*) 'gamma = ',gamma
      write(*,*) 'Energy = Quantised energy level En'
      write(*,*) '--------------------------------------------'
      write(*,*) '    step     tol        Energy         Error'
      write(*,*) '--------------------------------------------'

      do 100 i = 0, 5, +1
            result  = false_position(n, gamma, step, en, tol)
            error = dabs(answer - result)
            write(*,20) step, tol, result, error
20          format(i8, f9.7, f14.10, e14.6)
            step = step*10
            tol = 1.0d0/step
100   continue

      write(*,*) '--------------------------------------------'

      end

c-------------------------------------------------------------
c FUNCTION        : xin
c-------------------------------------------------------------

      double precision function xin(en)    
      implicit none
      real*8 en

      xin = (3.0d0-dsqrt(en + 1.0d0))/2.0d0
      return

      end

c-------------------------------------------------------------
c FUNCTION        : xout
c-------------------------------------------------------------

      double precision function xout(en)
      implicit none
      real*8 en

      xout = (3.0d0+dsqrt(en + 1.0d0))/2.0d0
      return

      end

c-------------------------------------------------------------
c FUNCTION        : f
c-------------------------------------------------------------

      double precision function f(en, n, gamma, step)
      implicit none
      integer*4 n, step
      real*8 xin, xout, simpson
      real*8 en, gamma, xi, xo, pi

      pi = 4*atan(1.0d0)
      xi = xin(en)
      xo = xout(en)

      f = gamma*simpson(step,xi,xo,en) - pi*(dfloat(n) + 0.5d0)
      return

      end

c-------------------------------------------------------------
c FUNCTION        : g
c PURPOSE         : integrand
c-------------------------------------------------------------

      double precision function g(en, x)
      implicit none
      real*8 en, x, v, arg

      arg = en - v(x)

      if (dabs(arg) .lt. 1.0E-14) then 
        g = 0.0d0 
      else 
        g = dsqrt(arg)
      endif
      return

      end

c-------------------------------------------------------------
c FUNCTION        : v   
c PURPOSE         : potential
c-------------------------------------------------------------


      double precision function v(x)
      implicit none
      real*8 x    

      v = 4*(x-1)*(x-2)
      return

      end

c-------------------------------------------------------------
c FUNCTION        : simpson
c-------------------------------------------------------------

      double precision function simpson(step, a, b, en)
      implicit none
      integer*4 step, i, factor
      real*8 a, b, en
      real*8 h, sum, x, g

      i = 0
      factor = 4
      sum = 0.0d0
      x = 0.0d0
      h = (b-a)/dfloat(step)

      do 300 i = 1, (step-1), +1
            x = a+i*h
            if (factor .eq. 2) then
                  sum = sum + 2.0d0*g(en,x)
                  factor = 4
            else
                  sum = sum + 4.0d0*g(en,x)
                  factor = 2
            endif
300   continue
      sum = sum + g(en,a) + g(en,b)
      sum = (h*sum)/3.0d0

      simpson = sum
      return

      end

c-------------------------------------------------------------
c FUNCTION  : false_position
c PURPOSE   : calculates and returns the root of the function
c           : f, after position x1, to a tolerance of tolx
c           : using the false position method. The parameter
c           : step determines the accuracy to which the function
c           : f can be calculated to.
c-------------------------------------------------------------

      double precision function false_position(n,gamma,step,
     > start,tolx)
      implicit none 
      integer*4 step, n 
      real*8 start, tolx, gamma 
      real*8 f 
      real*8 x1, x2, f1, f2, x3, f3, h 

      h = 0.3d0

      x1 = start 
      x2 = x1 + h 
      f1 = f(x1,n,gamma,step) 
      f2 = f(x2,n,gamma,step) 

      do while (f1*f2 .ge. 0.0d0) 
            x2 = x2 + h 
            f2 = f(x2,n,gamma,step) 
90          format(f15.8, f15.8)
      end do 
      
      x3 = x2 - f2*(x2-x1)/(f2-f1) 
      f3 = f(x3,n,gamma,step)

      do while (dabs(f3) .gt. tolx)
            if (f1*f3 .lt. 0.0d0) then
                  x2 = x3
            else
                  x1 = x3
            endif
            f1 = f(x1,n,gamma,step)
            f2 = f(x2,n,gamma,step)
            x3 = x2 - f2*(x2-x1)/(f2-f1)
            f3 = f(x3,n,gamma,step)
      end do

      false_position = x3
      return

      end

