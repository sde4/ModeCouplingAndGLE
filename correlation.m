function [ cf] = correlation(file1, file2, col)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tmp=importdata(file1);
tmp1 = tmp(:, col);
tmp=importdata(file2);
tmp2 = tmp(:, col);
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