C-------------------------------------------------------------
C PROGRAM      : main
C-------------------------------------------------------------

    program main 
    implicit none 
    integer i, n, step 
    real*8 gamma, en, answer, error, xmin, a
	real*8 result, tol, false_position

      write(*,*)'Enter n, a'
      read(*,*) n, a

      answer = -0.943121971d0
      gamma = 33.6567d0*a

      en = -1.0d0
      step = 10
      tol = 1.0d0/step
      xmin = 0.74166d0*a

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
            result  = false_position(n, gamma, step, en, tol, xmin)
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

      double precision function xin(en, xmin)    
      implicit none
      real*8 en, xmin

      xin = xmin - dlog(1.0d0 + dsqrt(en + 1.0d0))
      return

      end

c-------------------------------------------------------------
c FUNCTION        : xout
c-------------------------------------------------------------

      double precision function xout(en, xmin)
      implicit none
      real*8 en, xmin

c      write(*,*)'[xout]:[en] ',en
      xout = xmin - dlog(1.0d0 - dsqrt(en + 1.0d0))
      return

      end

c-------------------------------------------------------------
c FUNCTION        : f
c-------------------------------------------------------------

      double precision function f(en, n, gamma, step, xmin)
      implicit none
      integer*4 n, step
	  real*8 xin, xout, simpson
	  real*8 en, xmin, gamma, xi, xo, pi

      pi = 4*atan(1.0d0)
      xi = xin(en, xmin)
      xo = xout(en, xmin)

      f = gamma*simpson(step,xi,xo,en,xmin) - pi*(dfloat(n) + 0.5d0)
      return

      end

c-------------------------------------------------------------
c FUNCTION        : g
c PURPOSE         : integrand
c-------------------------------------------------------------

      double precision function g(en, x, xmin)
      implicit none
      real*8 en, xmin, x, v, arg

      arg = en - v(x, xmin)

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


      double precision function v(x, xmin)
      implicit none
      real*8 x, xmin

      v = (1 - dexp(-(x - xmin)))**2.0d0 - 1.0d0
      return

      end

c-------------------------------------------------------------
c FUNCTION        : simpson
c-------------------------------------------------------------

      double precision function simpson(step, a, b, en, xmin)
      implicit none
      integer*4 step, i, factor
      real*8 a, b, en, xmin
      real*8 h, sum, x, g

      i = 0
      factor = 4
      sum = 0.0d0
      x = 0.0d0
      h = (b-a)/dfloat(step)

      do 300 i = 1, (step-1), +1
            x = a+i*h
            if (factor .eq. 2) then
                  sum = sum + 2.0d0*g(en,x,xmin)
                  factor = 4
            else
                  sum = sum + 4.0d0*g(en,x,xmin)
                  factor = 2
            endif
300   continue
      sum = sum + g(en,a,xmin) + g(en,b,xmin)
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
     > start,tolx,xmin)
      implicit none 
      integer*4 step, n 
      real*8 start, tolx, gamma 
      real*8 f, xmin 
      real*8 x1, x2, f1, f2, x3, f3, h 

      h = 0.00001d0

      x1 = start 
      x2 = x1 + h 
      f1 = f(x1,n,gamma,step,xmin) 
      f2 = f(x2,n,gamma,step,xmin) 

      do while (f1*f2 .ge. 0.0d0) 
            x2 = x2 + h 
            f2 = f(x2,n,gamma,step,xmin) 
      end do 

      x3 = x2 - f2*(x2-x1)/(f2-f1) 
      f3 = f(x3,n,gamma,step,xmin)

      do while (dabs(f3) .gt. tolx)
            if (f1*f3 .lt. 0.0d0) then
                  x2 = x3
            else
                  x1 = x3
            endif
            f1 = f(x1,n,gamma,step,xmin)
            f2 = f(x2,n,gamma,step,xmin)
            x3 = x2 - f2*(x2-x1)/(f2-f1)
            f3 = f(x3,n,gamma,step,xmin)
      end do

      false_position = x3
      return

      end

