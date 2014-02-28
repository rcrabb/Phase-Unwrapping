%%

% P(B|D) and P(D|B) notes

% sample value of calibration constant I
I = 14.7113;
% Distance D is measured in meters
warning('off','MATLAB:quad:ImproperFcnValue');
% Range of values for B is [0 2.5]
Bs = [ .1 .2 .5 1 1.5 2.5 ];
figure('Color',[1.0 1.0 1.0]);
for iB = 1:length(Bs)
B = Bs(iB);

maxD = (I/B)^.5;

% distribution based on uniform prior for albedo and Beta
F = @(D)( D.^2./I .* log( B.*D.^2 ./ (I-(I.^2-B.^2*D.^4).^.5) ));
% get the numerical integration as close to the bounds as matlab will allow
delta = .0001;
intPb = NaN;
while (isnan(intPb))
    intPb = quad(F,delta,maxD-delta);
    delta = delta*10;
end

G = @(D) real( F(D)/intPb );


% distribution based on uniform prior for albedo and
% cos(beta)sin(beta) for surface normal Beta
H = @(D)( 2*D.^2./I .* (1 - B.*D.^2./I) );
% get the numerical integration as close to the bounds as matlab will allow
delta = .0001;
intPb = NaN;
while (isnan(intPb))
    intPb = quad(H,delta,maxD-delta);
    delta = delta*10;
end

J = @(D) ( H(D)/intPb );



% Check that function integrates to 1
delta = .0001;
intPb = NaN;
while (isnan(intPb))
    intPb = quad(G,delta,maxD-delta);
    delta = delta*10;
end
disp(['Integral for B = ' num2str(B) ' is ' num2str(intPb)]);

subplot(2,3,iB);fplot(@(x)[G(x),J(x)],[0 maxD 0 maxD]);xlim([0 15]);ylim([0 1]);
xlabel('D (m)');
%xlabel('D from 0 to \surd{ (Brightness / calib const) }');
ylabel('P(D|B)');
%ylabel('P(D|B) = P(B|D) / \int_D P(B|D)');
legend('p(\beta) \propto const','p(\beta) \propto sin(\beta)cos(\beta)');

sTitle = ['P(D|B) for Brightness val ' num2str(B)];
title(sTitle);
end
%% Let's try for P(D|B) specifically

% P(B|D) and P(D|B) notes

% sample value of calibration constant I
I = 14.7113;
% Distance D is measured in meters

% Range of values for B is [0 2.5]
Bs = [ .1 .2 .5 1 1.5 2.5 ];
Dlim = (I/Bs(1))^.5;
figure('Color',[1.0 1.0 1.0]);
for iB = 1:length(Bs)
B = Bs(iB);

maxD = (I/B)^.5;

F = @(D)( 2*D.*B./I .* (1 - B.*D.^2./I) );

% get the numerical integration as close to the bounds as matlab will allow
delta = .0001;
intPb = NaN;
while (isnan(intPb))
    intPb = quad(F,delta,maxD-delta);
    delta = delta*10;
end

disp(['Integral for P(D|B) = ' num2str(B) ' is ' num2str(intPb)]);
G = @(D) ( F(D)/intPb );

% Check that function integrates to 1
delta = .0001;
intPb = NaN;
while (isnan(intPb))
    intPb = quad(G,delta,maxD-delta);
    delta = delta*10;
end
disp(['Integral for B = ' num2str(B) ' is ' num2str(intPb)]);

subplot(2,3,iB);fplot(G,[0 maxD]);xlim([0 ceil(Dlim)]);ylim([0 0.7]);
xlabel('D (m)');
%xlabel('D from 0 to \surd{ (Brightness / calib const) }');
ylabel('P(D|B) = 4BD/I (1 - BD^2/I');

sTitle = ['P(D|B) for Brightness val ' num2str(B)];
title(sTitle);

% capture the data to be plotted for a group plot
input{iB} = 0:.01:maxD;
output{iB} = G(input{iB});
end

figure('Color',[1.0 1.0 1.0]);
plot(input{1},output{1},'-',input{3},output{3},'--',input{5},output{5},':',input{6},output{6},'-.');
legend('B = 0.1','B = 0.5','B = 1.5','B = 2.5');
xlabel('D (m)');
ylabel('P(D|B)')
title('P(D|B) for various Brightness values');

%% Let's try for P(B|D) specifically

% P(B|D) and P(D|B) notes

% sample value of calibration constant I
I = 14.7113;
% Distance D is measured in meters

% Range of values for D is [0 8]
Ds = [ 1 2 3 4 5 6];
figure('Color',[1.0 1.0 1.0]);
for iD = 1:length(Ds)
D = Ds(iD);

maxB = (I./D.^2)
maxP = 2*D.^2./I

F = @(B)( 2*D.^2./I .* (1 - B.*D.^2./I) );

% get the numerical integration as close to the bounds as matlab will allow
delta = 1e-10;
intPb = NaN;
while (isnan(intPb))
    intPb = quad(F,delta,maxB-delta);
    delta = delta*10;
end
delta
intPb
G = @(D) ( F(D)/intPb );

% Check that function integrates to 1
delta = 1e-10;
intPb = NaN;
while (isnan(intPb))
    intPb = quad(G,delta,maxB-delta);
    delta = delta*10;
end
disp(['Integral for D = ' num2str(D) ' is ' num2str(intPb)]);

subplot(2,3,iD);fplot(G,[0 maxB]);xlim([0 15]);ylim([0 2]);
xlabel('Brightness');
%xlabel('D from 0 to \surd{ (Brightness / calib const) }');
ylabel('P(B|D)');

sTitle = ['P(B|D) for Distance ' num2str(D) 'm'];
title(sTitle);
end
%% The path integral loop!
iter = 0;
for c = .00005:.0001:1
iter = iter+1;
F = @(x)(1+(sec(x).*tan(x).*c).^2).^.5;
Q(iter) = quad(F,0,acos(c));
end

%%

%% The path integral loop!
iter = 0;
%for c = .00005:.0001:1
for c = .0005:.001:1
iter = iter+1;
F = @(x)(sin(acos(c./x))./x);
P(iter) = c*quad(F,c,1);
end
figure;plot(.0005:.001:1,P)
xlabel('c');
ylabel('\int_c^1 (sin(cos^{-1}(c/x))/x)dx')

%% The path integral loop!
iter = 0;
%for c = .00005:.0001:1
for c = .0005:.001:1
iter = iter+1;
F = @(x)(sin(x).*cos(x));
O(iter) = quad(F,0,acos(c));
end
figure;plot(.0005:.001:1,O)
xlabel('c');
%ylabel('$\frac{BD^2}{I}\int_0^{acos(c)}(sin(\beta)cos(\beta)d\beta$')
ylabel('c \int_0^{acos(c)}(sin(\beta)cos(\beta)d\beta')



%% Solved without path integral
iter = 0;
%for c = .00005:.0001:1
for c = .0005:.001:1
iter = iter+1;
F = @(x)( x.*(log(x)-log((1-x^2).^-.5-1)) );
R(iter) = F(c);
end
figure;plot(.0005:.001:1,R)
xlabel('c');
%ylabel('$\frac{BD^2}{I}\int_0^{acos(c)}(sin(\beta)cos(\beta)d\beta$')
ylabel('c x.*(log(x)-log((1-x^2).^-.5-1))')

