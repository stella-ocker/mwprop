      subroutine getarg(arg,argt)

c     get a command line argument using gfortran

      implicit none

c input

c     number of argument to get
      integer arg

c output

c     string
      character*(*) argt

c executable

      call get_command_argument(arg,argt)

      return
      end
