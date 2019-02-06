
using Base.Cartesian
using CuArrays
for N = 1:5
    @eval begin
        function laplacian(A::Array{T,$N}) where T<:Number
            B = similar(A)
            @nloops $N i A begin
                tmp = zero(T)
                # Do the shift by +1.
                @nexprs $N d1->begin
                    tmp += (i_d1 < size(A,d1)) ? (@nref $N A d2->(d2==d1) ? i_d2+1 : i_d2) : (@nref $N A i)
                end
                # Do the shift by -1.
                @nexprs $N d1->begin
                    tmp += (i_d1 > 1) ? (@nref $N A d2->(d2==d1) ? i_d2-1 : i_d2) : (@nref $N A i)
                end
                # Subtract the center and store the result
                (@nref $N B i) = tmp - 2*$N*(@nref $N A i)
            end
            B
        end
    end
end






function δ⁻(size::Int, step::Number)
    δ = zeros(size,size)
    for i in 1:size
        δ[i,i] = 1
        if i > 1
            δ[i-1,i] = -1
        end
    end
    return (1/step).*δ
end


function δ⁺(size::Int, step::Number)
    δ = zeros(size,size)
    for i in 1:size
        δ[i,i] = -1
        if i < size
            δ[i+1,i] = 1
        end
    end
    return (1/step).*δ
end

function δ(size, step)
    δ = zeros(size,size)
    for i in 1:size
        if i > 1
            δ[i-1,i] = -1
        end
        if i < size
            δ[i+1,i] = 1
        end
    end
    return δ
end

function δδ(size::Int, step::Number)
     δ = δ⁻(size, step)*δ⁺(size, step)
     δ[2,1] = δ[2,1]*2
     δ[end-1,end] = δ[end-1,end]*2
     δ[end,end]=δ[1,1]
     return δ
  end

 function ⊗(A::Array{N,2},B::Array{N,2}, direction::Int) where N<:Number
     s = size(B)
     C = zeros(N,s)
     if direction == 1
         for i in 1:s[1]
             C[:,i] = A*B[:,i]
         end
    elseif direction == 2
        Bt = Transpose(B)
        C = Transpose(C)
        for i in 1:s[2]
            C[:,i] = A*Bt[:,i]
        end
        C = Transpose(C)
    else
        throw("direction too large, choose 1 for through columns and 2 for through rows")
    end
    return C
end

function ⊗(A::Array{N,2},B::Array{Complex{Float64},2}, direction::Int) where N<:Number
    s = size(B)
    C = Complex.(zeros(s))
    if direction == 1
        for i in 1:s[1]
            C[:,i] = A*B[:,i]
        end
   elseif direction == 2
       Bt = Transpose(B)
       C = Transpose(C)
       for i in 1:s[2]
           C[:,i] = A*Bt[:,i]
       end
       C = Transpose(C)
   else
       throw("direction too large, choose 1 for through columns and 2 for through rows")
   end
   return C
end

function ⊗(A::CuArray,B::CuArray, direction::Int64)
    s = size(B)
    C = cu(zeros(s))
    if direction == 1
        for i in 1:s[2]
            @inbounds C[:,i] = A*B[:,i]
        end
   elseif direction == 2
       for i in 1:s[1]
           @inbounds C[i,:] = A*B[i,:]
       end
   end
   return C
end

∇(u, δxx, δyy) = ⊗(δxx,u,1).+⊗(δyy,u,2)
