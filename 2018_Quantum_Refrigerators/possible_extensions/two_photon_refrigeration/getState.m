ket00 = kron([1 0 0]',[1 0 0]');
ket01 = kron([1 0 0]',[0 1 0]');
ket02 = kron([1 0 0]',[0 0 1]');
ket10 = kron([0 1 0]',[1 0 0]');
ket11 = kron([0 1 0]',[0 1 0]');
ket12 = kron([0 1 0]',[0 0 1]');
ket20 = kron([0 0 1]',[1 0 0]');
ket21 = kron([0 0 1]',[0 1 0]');
ket22 = kron([0 0 1]',[0 0 1]');
rho1 = ket00*(ket00+ket02+ket20+ket22+ket11)' + ket02*(ket00+ket02+ket20+ket22+ket11)' + ...
 ket20*(ket00+ket02+ket20+ket22+ket11)' +ket22*(ket00+ket02+ket20+ket22+ket11)' +...
ket11*(ket00+ket02+ket20+ket22+ket11)';

rho2 = ket01*(ket01+ket10+ket21+ket12)' + ket10*(ket01+ket10+ket21+ket12)' + ...
ket12*(ket01+ket10+ket21+ket12)'+ket21*(ket01+ket10+ket21+ket12)';
