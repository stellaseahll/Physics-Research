# some helper functions

#save a 2D-array in comma-separated form
function saveMatrix(M, filename::String, header::String = "")
  f = open(filename, "w")
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
