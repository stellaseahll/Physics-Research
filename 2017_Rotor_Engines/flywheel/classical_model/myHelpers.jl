# some helper functions

#save a 2D-array in comma-separated form
function saveMatrix(M, filename::String, header::String = "")
  f = open(filename, "w")
  print(f, "# created on ")
  println(f, string(now()))
  println(f, header)
  (z,s) = size(M)
  for i=1:z
    for j=1:s
      print(f, M[i,j])
      if j<s
        print(f, ", ")
      else
        println(f, "")
      end
    end
  end
  close(f)
  nothing
end

##  test routines

#Compare cutoff for negative numbers with if or max
# --> max seems faster!
function randNegCut!(A)
  @inbounds @simd for i=1:length(A)
    A[i] = randn()
    if A[i]<0
      A[i]=0.0
    end
  end
end
function randNegCut2!(A)
  @inbounds @simd for i=1:length(A)
    A[i] = max(0.0,randn())
  end
end
