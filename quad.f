c-------------------------------------------------------------
c PROGRAM   : QUADRATURE - TRAPEZOIDAL AND SIMPSON'S RULES
c 
c SUBJECT   : 640-364 COMPUTATIONAL PHYSICS
c NAME      : MICHAEL PAPASIMEON
c DATE      : 30/07/1997
c
c PURPOSE   : This program does a numerical integration 
c           : between 0 and 1 for the function ln(1+x)/x
c           : using both the Trapezoidal and Simpson's
c           : rules for quadrature.
c-------------------------------------------------------------
	
      program main
            call all_trapezoidals
            call all_simpsons
      end

c-------------------------------------------------------------
c FUNCTION  : f
c PURPOSE   : Given a real number x, returns the result of
c           : evaluating ln(1+x)/x.
c           : Care is taken when x = 0, since this causes
c           : a division by zero error, and hence a value
c           : of 1 is returned in this case.
c-------------------------------------------------------------

      double precision function f(x)
      implicit none
      real*8 x
      if ( x .eq. 0.0d0 ) then
            f = 1.0
      else
            f = dlog(1+x)/x
      endif
      return 
      end

c-------------------------------------------------------------
c SUBROUTINE : all_trapezoidals
c PURPOSE    : Calls the trapezoidal subroutine 6 each time
c            : with an increase in the number of divisions
c            : (increased by a factor of 10). The result
c            : is a table of the result of the quadrature
c            : including an output of the percentage error.
c-------------------------------------------------------------

      subroutine all_trapezoidals 
      implicit none 
      integer*4 n, i
      real*8 actual, pi, a, b

      a = 0.0d0
      b = 1.0d0
      pi = 4*atan(1.0d0)
      actual = (pi**2)/dfloat(12)

      write(*,*) '---------------------------------------------------'
      write(*,*) '                TRAPEZOIDAL RULE'
      write(*,*) '---------------------------------------------------'
      write(*,*) 'Lower Limit : a = ',a
      write(*,*) 'Upper Limit : b = ',b
      write(*,*) 'Actual Result = ',actual
      write(*,*) '---------------------------------------------------'
      write(*,*) '        n             h        Result         Error'
      write(*,*) '---------------------------------------------------'

      n = 10
      do 10 i = 0, 5, +1
            call trapezoidal(n,a,b,actual) 
            n = n*10
10    continue
      write(*,*) '---------------------------------------------------'

      end 
     
c-------------------------------------------------------------
c SUBROUTINE : trapezoidal
c PURPOSE    : Performs numerical quadrature using the 
c            : trapezoidal rule for the function f
c            : between a and b and prints out the result
c            : for the number of divisions given as a 
c            : parameter.
c-------------------------------------------------------------

      subroutine trapezoidal(n, a, b, answer)
      implicit none 
      integer*4 n
      real*8 a, b, answer
	integer*4 i
      real*8 error, h, sum, x, f 

	i = 0
      sum = 0.0d0
      x = 0.0d0
      h = (b-a)/dfloat(n)

	do 100 i = 1, (n-1), +1
		x = a + i*h
		sum = sum + 2.0d0*f(x)
100   continue
	sum = sum + f(a) + f(b)
      sum = (0.5d0)*h*sum
      error = abs(sum - answer)

      write(*,200) n, h, sum, error
200   format(i10, f14.10, f14.10, g14.6)

      end

c-------------------------------------------------------------
c SUBROUTINE : all_simpsons
c PURPOSE    : Calls the simpson subroutine 6 each time
c            : with an increase in the number of divisions
c            : (increased by a factor of 10). The result
c            : is a table of the result of the quadrature
c            : including an output of the percentage error.
c-------------------------------------------------------------

      subroutine all_simpsons
      implicit none 
      integer*4 n, i
      real*8 actual, pi, a, b

      a = 0.0d0
      b = 1.0d0
      pi = 4*atan(1.0d0)
      actual = (pi**2)/dfloat(12)

      write(*,*) '---------------------------------------------------'
      write(*,*) '                SIMPSONS RULE'
      write(*,*) '---------------------------------------------------'
      write(*,*) 'Lower Limit : a = ',a
      write(*,*) 'Upper Limit : b = ',b
      write(*,*) 'Actual Result = ',actual
      write(*,*) '---------------------------------------------------'
      write(*,*) '        n             h        Result         Error'
      write(*,*) '---------------------------------------------------'

      n = 10 
      do 10 i = 0, 5, +1
            call simpson(n,a,b,actual) 
            n = n*10
10    continue
      write(*,*) '---------------------------------------------------'

      end 

c-------------------------------------------------------------
c SUBROUTINE : simpson
c PURPOSE    :
c-------------------------------------------------------------

      subroutine simpson(n, a, b, answer)
      implicit none
      integer*4 n, i, factor
      real*8 a, b, answer
      real*8 error, h, sum, x, f

      i = 0
      factor = 4
      sum = 0.0d0
      x = 0.0d0
      h = (b-a)/dfloat(n)

	do 300 i = 1, (n-1), +1
		x = a+i*h
		if (factor .eq. 2) then 
			sum = sum + 2*f(x) 
			factor = 4 
		else
            	sum = sum + 4*f(x) 
			factor = 2
		endif
300   continue
      sum = sum + f(a) + f(b)
      sum = (h*sum)/3.0d0
      error = dabs(sum - answer)

      write(*,400) n, h, sum, error
400   format(i10, f14.10, f14.10, g14.6) 

      end
      

