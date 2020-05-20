function [H, BW] = channelParameters(time, h_t)

H = 0;
M1 = 0;
M2 = 0;

for I = 1:length(time)
   H = H + h_t(I);
   M1 = M1 + h_t(I)*time(I);
   M2 = M2 + h_t(I)*time(I)^2;
end

M2 = M2/H;
M1 = M1/H;

BW = 0.2/sqrt(M2 - M1^2);