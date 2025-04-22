struct FieldMapTracking end
function _track!(
  i,
  v,
  work,
  bunch::Bunch,
  ele::Union{LineElement,BitsLineElement}, 
  ::FieldMapTracking,
)
  # Unpack the line element
  ma = ele.AlignmentParams
  bm = ele.BMultipoleParams
  bp = ele.BendParams
  L = ele.L

  # Function barrier
  linear_universal!(i, v, work, bunch, L, bm, bp, ma)
end