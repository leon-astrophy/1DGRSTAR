F90=pgfortran
F90FLAGS= -fastsse -fast -Mipa=fast,inline -Msmartalloc -Mfprelaxed -Mstack_arrays -O4 -Mvect=prefetch -Mprefetch -mcmodel=medium -Minfo
MODINC="-I ./"
LDFLAGS= -fastsse -fast -Mipa=fast,inline -Msmartalloc -Mfprelaxed -Mstack_arrays -O4 -Mvect=prefetch -Mprefetch -mcmodel=medium -Minfo