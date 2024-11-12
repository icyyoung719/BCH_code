% generateSi_a takes in a syndrome component function, S_i(x), and a power
% of alpha, i, and returns S_i(alpha^i).
function [ Si_a ] = generateSi_a(Si_x, i, m)
    length_of_Si_a = 2^m;
	Si_a = de2bi(0,length_of_Si_a);
	temp = bi2de(Si_x, 'left-msb');
	Si_x = de2bi(temp);

	% Multiply all bits, except for the lsb, by i.
	for j = 2:length(Si_x)
		Si_x(j) = Si_x(j)*i;
	end

	for k = 2:length(Si_x)
		Si_x(k) = Si_x(k)*(k-1);
	end

	for k = 2:length(Si_x)
		if Si_x(k) ~= 0  
			Si_a(length_of_Si_a-Si_x(k)) = 1; 
		end
	end

	if Si_x(1) == 1
		Si_a(length_of_Si_a) = 1;
end
