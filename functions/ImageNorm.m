function Output=ImageNorm(Input)
Input=double(Input);
MinInput=min(Input(:));
MaxInput=max(Input(:));
InvertRange=1/(MaxInput-MinInput+eps);
Output=double((Input-MinInput)*InvertRange);
