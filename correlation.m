function [ cf] = correlation(file1, file2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tmp1=importdata(file1); 
tmp2=importdata(file2);
len1 = length(tmp1);
len2 = length(tmp2);
if len1<len2
  len=len1;
  tmp2(len+1:end)=[];
elseif len2<len1
  len=len2;
  tmp1(len+1:end)=[];
else
  len=len1;
end
tmp      = xcorr(tmp1-mean(tmp1), tmp2-mean(tmp2),  'unbiased');
cf = tmp(len:len+round(len/10));

end

