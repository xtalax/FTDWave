module PDEUtils
export  δ⁻, δδ, ∇

using Base.Cartesian
using LinearAlgebra
using DiscreteAxis
using SparseArrays
#using CuArrays
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

#function laplacian(A::AbstractArray{T}, ) where T <: Number





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

function δ⁺(F::AbstractVector{N}, axis::DAxis) where N <: Number
    δ = similar(F)
    for i in axis.i[1:end-1]
        δ[i] = (F[i+1]-F[i])/axis.Δ[i]
    end
    δ[end] = -F[end]
    return δ
end

function δ⁻(F::AbstractVector{N}, axis::DAxis) where N <: Number
    δ = similar(F)
    for i in axis.i[2:end]
        δ[i] = (F[i]-F[i-1])/axis.Δ[i]
    end
    return δ
end

function δδ(F::AbstractVector{N}, axis::DAxis) where N <: Number
    δ = similar(F)
    for i in axis.i[2:end-1]
        δ[i] = 2*F[i]-F[i-1]-F[i-2]/axis
    end
    return δ
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
    return sparse(δ)
end

function δδ(size::Int, step::Number)
     δ = zeros(size,size)
     for i in 1:size
         δ[i,i] = -2.0

         if i >= 2
             δ[i-1,i] = δ[i,i-1] = 1.0
         end
     end
     δ[1,1:4] =  [2.0,-5.0,4.0,-1.0]
      δ[end, end-3:end] = [-1.0,4.0,-5.0,2.0]
     return δ./step^2
  end

 function ⊗(A::AbstractArray{N,2},B::AbstractArray{N,2}, direction::Int) where N<:Number
     s = size(B)
     C = zeros(N,s)
     if direction == 1
         for i in 1:s[1]
             @views C[:,i] = A*B[:,i]
         end
    elseif direction == 2
        Bt = Transpose(B)

        for i in 1:s[2]
            @views C[:,i] = A*Bt[:,i]
        end
        C = Transpose(C)
    else
        throw("direction too large, choose 1 for through columns and 2 for through rows")
    end
    return C
end

 function SliceProduct(A::Array{N,2},B::Array{N,2}, direction::Int) where N<:Number
     s = size(B)
     C = zeros(N,s)
     if direction == 1
         for i in 1:s[1]
            @views C[:,i] = A*B[:,i]
         end
    elseif direction == 2
        for i in 1:s[2]
            @views C[i,:] = A*B[i,:]
        end
    else
        throw("direction too large, choose 1 for through columns and 2 for through rows")
    end
    return C
end
#=
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
=#
∇(u, δxx, δyy) = ⊗(δxx,u,1).+⊗(δyy,u,2)
end
