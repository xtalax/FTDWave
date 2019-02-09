module SignalProcessing
export rect, convolve, convolvepower, unitize, correlation
using FFTW
FFTW.set_num_threads(4)

rect(x::Number) = abs(x) > 0.5 ? 0 : 1

function convolve(a::AbstractArray,b=a::AbstractArray)
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

function correlation(a,b=a)
    convolve(a,b[end:-1:1])
end

unitize(A::AbstractArray) = A./(maximum(abs.(A)))

dB(a::Number) = 20*log10(abs(a))

#x = -2:0.01:2

#fs = 1/0.01
#fgrid = fs*(0:(length(x)-1))/(length(x))

#Ft = rect.(x)
#F = rfft(Ft)/length(x)
#plot(x,Ft)
#plot(fgrid,F)
end
