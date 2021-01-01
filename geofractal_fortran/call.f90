program call_geofractal
use types; 
implicit none
!--------------------------------------------------------------------------------
!       some cluster constants 
!--------------------------------------------------------------------------------
real(kind=dp),parameter :: k0_chain = sqrt(3.0_dp)
real(kind=dp),parameter :: df_chain = 1.0_dp
real(kind=dp),parameter :: df_bpca  = 3.0_dp
real(kind=dp),parameter :: k0_bpca  = 0.3_dp
real(kind=dp),parameter :: k0_bcca  = 1.04_dp
real(kind=dp),parameter :: df_bcca  = 1.90_dp

!--------------------------------------------------------------------------------

integer,parameter       :: n        = 250
integer                 :: i,iqapp,iqcon,iqcor
real(kind=dp)           :: PN(1:n),G(1:n),G0,k0,df

!--------------------------------------------------------------------------------
!       calculation mode 
!--------------------------------------------------------------------------------
iqapp   = 3                     ! switch 1: 1 or 3
iqcon   = 2                     ! switch 2: 1 or 2
iqcor   = 3                     ! correlation func.
!--------------------------------------------------------------------------------
!       straight chain 
!--------------------------------------------------------------------------------
!k0      = k0_chain              ! fractal prefactor
!df      = df_chain              ! fractal dimension
!--------------------------------------------------------------------------------
!       BPCA
!--------------------------------------------------------------------------------
!df      = df_bpca               ! fractal dimension
!k0      = k0_bpca               ! fractal prefactor
!--------------------------------------------------------------------------------
!       BCCA
!--------------------------------------------------------------------------------
!df      = df_bcca               ! fractal dimension
!k0      = k0_bcca               ! fractal prefactor
!--------------------------------------------------------------------------------
!       other fractal dimensions
!--------------------------------------------------------------------------------
df      = 3.0_dp                ! fractal dimension
k0      = (k0_chain-k0_bpca)/(df_chain-df_bpca)*(df-df_bpca)+k0_bpca
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
!       call geofractal
!--------------------------------------------------------------------------------
do i=1,n
        PN(i) = 1.0_dp * (1.e10_dp/1.0_dp) ** &
                (real(i-1,kind=dp)/real(n-1,kind=dp))
        call geofrac(iqapp,iqcon,iqcor,PN(i),k0,df,G0)
        G(i) = G0
enddo

!--------------------------------------------------------------------------------
!       Output
!--------------------------------------------------------------------------------
open(10,file="gratio.out",status="unknown")
write(10,1100) iqapp," = iqapp"
write(10,1100) iqcon," = iqcon"
write(10,1100) iqcor," = iqcor"
write(10,1000) df, " = fractal dimension"
write(10,1000) k0, " = fractal prefactor"
write(10,2000) "PN","G/NpiR0^2"
do i=1,n
        write(10,3000) PN(i),G(i)
enddo
stop

1000 format('#',1PE15.5,A)
1100 format('#',I15,A)
2000 format('#',8A15)
3000 format(' ',1P8E15.5)
print *, 'finished.'

stop
end program call_geofractal
