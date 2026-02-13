      function iargc()

c     get number of command line arguments using gfortran

      implicit none

c output

c     number of command line arguments
      integer iargc

c executable

      iargc=command_argument_count()

      return
      end
