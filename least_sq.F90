!This file is part of TDEFit.

!TDEFit is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TDEFit is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with TDEFit.  If not, see <http://www.gnu.org/licenses/>.

SUBROUTINE least_sq(x,y,a,b,siga,sigb,chi2,q,sig)
    USE tdefit_util, ONLY : gamma_inc, assert_eq
    real, DIMENSION(:), INTENT(IN) :: x,y
    real, INTENT(OUT) :: a,b,siga,sigb,chi2,q
    real, DIMENSION(:), OPTIONAL, INTENT(IN) :: sig
    integer :: ndata
    real :: sigdat,ss,sx,sxoss,sy,st2
    real, DIMENSION(size(x)), TARGET :: t
    real, DIMENSION(:), POINTER :: wt
    if (present(sig)) then
        ndata=assert_eq(size(x),size(y),size(sig),'least squares')
        wt=>t
        wt(:)=1.0/(sig(:)**2)
        ss=sum(wt(:))
        sx=dot_product(wt,x)
        sy=dot_product(wt,y)
    else
        ndata=assert_eq(size(x),size(y),'fit')
        ss=real(size(x))
        sx=sum(x)
        sy=sum(y)
    end if
    sxoss=sx/ss
    t(:)=x(:)-sxoss
    if (present(sig)) then
        t(:)=t(:)/sig(:)
        b=dot_product(t/sig,y)
    else
        b=dot_product(t,y)
    end if
    st2=dot_product(t,t)
    b=b/st2
    a=(sy-sx*b)/ss
    siga=sqrt((1.0+sx*sx/(ss*st2))/ss)
    sigb=sqrt(1.0/st2)
    t(:)=y(:)-a-b*x(:)
    q=1.0
    if (present(sig)) then
        t(:)=t(:)/sig(:)
        chi2=dot_product(t,t)
        if (ndata > 2) q=gamma_inc(0.5*(size(x)-2),0.5*chi2)
    else
        chi2=dot_product(t,t)
        sigdat=sqrt(chi2/(size(x)-2))
        siga=siga*sigdat
        sigb=sigb*sigdat
    end if
END SUBROUTINE least_sq
