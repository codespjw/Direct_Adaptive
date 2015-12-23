m1 = IMG.vcm;
m2 = IMG.ma;
m2 = [ones(1,115) zeros(1,173-115); m2];

m = m1;
n = size(m2,1);
m(end+1:end+n,:)= [repmat(m1(end,1:58),n,1) m2(:,59:173)];
m3 = m(:,58)==m(:,59);
for i =1:size(m3,1)
    if m3(i) == 1
        m(i,59:end)=~m(i,59:end);
    end
end

i =1;
% while(i<size(m,1))
%     if all(m(i,:)==m(i+1,:))
%         m(i+1,:)=[];
%     else
%         i=i+1;
%     end
% end
m = [ones(1,58) zeros(1,173-58); m];
b = m(:,1:end-1)~=m(:,2:end);
b = double(b);
shade = 1.1;
c = b;
for j = 2:size(b,1)
    if all(c(j,:)==c(j-1,:))
        f = e;
        e = find(b(j,f+1:end),1)+f;
        if isempty(e)
            e = period/2-1;
        end
    else
        e = find(b(j-1,:)~=1 & b(j,:)==1,1);
        f = find(b(j,1:e-1));
    end
    if isempty(f)
        f = 1;
    end
    f = f(end);
    b(j,f+(f~=1):e-1) = shade;
end
b(1,1:57) = shade;
b(3:end+1,:) = b(2:end,:);
b(3,b(3,:)==shade)=0;
b(3,30:57)=shade;
% b(b==1) = 0.5;
b(b==0) = 2.3;
b(end,:)=[];
b(:,58) = 1.5;
figure; 
imagesc(b)

hold all
line(repmat([1,period/2-1],size(m,1),1)',0.5+[1:size(m,1);1:size(m,1)],'color','r','linewidth',2)
% colormap(hot)
ylabel('Iteration');
xlabel('Harmonic');
saveImgPdf(6,4,'frequency-paritions',1);