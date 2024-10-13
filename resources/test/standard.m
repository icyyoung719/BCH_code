%BCH码字长度
m=4;
n=2^m-1;
%消息长度
k=5;
%消息比特的行数
N=100;
%消息比特
msg=randi([0 1],N, k);
%BCH码的生成多项式
[genpoly,t]=bchgenpoly(n, k);
%BCH编码
code=bchenc(gf(msg), n,k);
%码字加入不超过纠错能力的误码
noisycode=code+randerr(N, n,1:t);
%BCH译码
[newmsg, err, ccode]=bchdec(noisycode, n, k);
if ccode==code
    disp('所有错误比特都被纠正')
end
if newmsg==msg
    disp('译码消息与原消息相同')
end
