
# ------------------------------------------------
# Codes I am using in Geoph 426-526
# M D Sacchi
# ------------------------------------------------

function SeisConvMatrix(w,n)
    # Make Toeplitz convolution Matrix such data y = conv (x,w) = W.x
    # Matlab has something call convmtx that does this.. 
    m = length(w);
    W = zeros(m+n-1,n)
    for k = 1:n
        W[k:k+m-1,k]=w
    end
    return W
end


function convolution(x,h)
    # convolve x with h
    # in Matlab this equivalent to "conv"

nx = length(x)
nh = length(h)
 y = zeros(nx+nh-1)

for ix = 1:nx
    for ih=1:nh
        y[ix+ih-1] = y[ix+ih-1] + x[ix]*h[ih]
    end
end

return y
end



function ls_decon(w,nf,mu)
    # Least-squares inverse filter of the wavelet w
    # nf is the filter lenght and mu thetrade-off 

  nw = length(w)

# Desired output 
   o = zeros(nw+nf-1); o[1] = 1. 

# set convolution matrix, right hand side term 
# and autocorrelation matrix R
 W = SeisConvMatrix(w,nf)
 R = W'*W;
 g = W'*o
 I = diagm(0=>ones(nf))

# least-squares solution 
 f = (R+mu*I)\g

# Actual output 
    oh = W*f

# Error 
    e = oh-o

return f,oh,e
end


function spectrum(s,dt) 
   # Compute spectrum of the signal s 
   # dt is sampling rate

 N = length(s);
 M = 501

# Freq. in radians in [-pi,pi]

 omega = [ -pi + 2.0*pi*(k-1)/(M-1) for k=1:M]
 f = omega/(2*pi*dt) 


 S =zeros(Complex,M); 
for k = 1: M 
   for  n = 1:N
    S[k] = S[k] + s[n]*exp(-im*omega[k]*(n-1))
   end
end

return omega, f, abs.(S)
end
