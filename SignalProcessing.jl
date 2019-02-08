module SignalProcessing
export rect, convolve, convolvepower, unitize
using FFTW
FFTW.set_num_threads(4)

rect(x::Number) = abs(x) > 0.5 ? 0 : 1

function convolve(a::AbstractVector,b=a::AbstractVector)
     X = fft(a)/length(a);  Y = fft(b)/length(b)
     Z = X.*Y
     return length(Z)*ifft(Z), floor(length(X)/2)
 end

function convolvepower(in,n)
    X = fft(in)
    X = X.^n

    Y = ifft(X)

    return Y
end
unitize(A::AbstractVecOrMat) = A./(maximum(abs.(A)))

#x = -2:0.01:2

#fs = 1/0.01
#fgrid = fs*(0:(length(x)-1))/(length(x))

#Ft = rect.(x)
#F = rfft(Ft)/length(x)
#plot(x,Ft)
#plot(fgrid,F)
end
