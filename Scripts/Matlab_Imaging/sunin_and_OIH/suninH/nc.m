function B=nc(A);

% combination of 'norm_to_uint8' and 'OIClip'

B=norm_to_uint8(OIClip(A, 1,4));
return;