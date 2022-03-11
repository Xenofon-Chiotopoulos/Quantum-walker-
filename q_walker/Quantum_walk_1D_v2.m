
%David Faux January 2021
%
% 1D line: set up for a "Fauxian" operator
%

globalForQuantumWalk1D_v2

fid_qubit=fopen(SummaryOutput_qubit,'wt');       % open output with identifier complex(a,b)
fid_probs=fopen(SummaryOutput_probs,'wt');       % open output with identifier pdf
fid_meanZ=fopen(SummaryOutput_meanZ,'wt');       % open output with identifier - ave(z posn)

F;
z = zeros(ntimes,linelength);
z2 = zeros(ntimes,linelength);
time = zeros(ntimes,1);
wt = zeros(ntimes,1);
vel_mn = zeros(ntimes,1);
vel2_mn = zeros(ntimes,1);
stdev_v = zeros(ntimes,1);
mean_z = zeros(ntimes,1);
mean_z2 = zeros(ntimes,1);
stdev_z = zeros(ntimes,1);

x = zeros(2,1);
y = zeros(2,1);
first = zeros(2,1);
second = zeros(2,1);

counter = zeros(ntimes,1);
detection_point = zeros(1,linelength);
detection_point(1,detector_point) = 1;

qur = zeros(ntimes,linelength,2);               % qubit array

qu = complex(qur);                              % the complex qur 
qu_y = complex(qur);
qu_first = complex(qur);
qu_second = complex(qur);

qunew = zeros(ntimes,linelength,2);             % array of new qubits at next time step
qunew_y = zeros(ntimes,linelength,2);
qunew_first = zeros(ntimes,linelength,2);
qunew_second = zeros(ntimes,linelength,2);

zprob = zeros(linelength,1);                    % product of prob*z
z2prob = zeros(linelength,1);                   % product of prob*z2
vel_prob = zeros(linelength,1);
vel2_prob = zeros(linelength,1);

prob = zeros(ntimes,linelength);
prob_y = zeros(ntimes,linelength);
prob_first = zeros(ntimes,linelength);
prob_second = zeros(ntimes,linelength);
prob_sum = zeros(ntimes,1);
probLine = zeros(linelength,1);

vel_mn = zeros(ntimes,1);
vel2_mn = zeros(ntimes,1);
sum_x = zeros(ntimes,1);
sum_x2 = zeros(ntimes,1);
uncertainty = zeros(ntimes,1);

Universe = zeros(ntimes,linelength);

% z coordinate with center as 0
for i = 1:linelength
    z(:,i) = i-linecenter;
end

% H = (sqrt(0.5))*[exp(1i*theta), 1; 1, exp(1i*theta)];           % the operator

%Hq = H*init_qubit;

qu(1,linecenter,1) = init_qubit(1,1);                            % creates the qu array at t=0
qu(1,linecenter,2) = init_qubit(2,1);                             % creates the qu array at t=0
qu_y(1,linecenter,1) = init_qubit_y(1,1);
qu_y(1,linecenter,2) = init_qubit_y(2,1);
qu_first(1,linecenter,1) = init_qubit_first(1,1);
qu_first(1,linecenter,2) = init_qubit_first(2,1);
qu_second(1,linecenter,1) = init_qubit_second(1,1);
qu_second(1,linecenter,2) = init_qubit_second(2,1);

prob(1,linecenter) = 1;                                         % prob at linecenter at t=0
prob_y(1,linecenter) = 1;
prob_first(1,linecenter) = 1;
prob_second(1,linecenter) = 1;

for t = 2:ntimes                                % t=0 is actually t=1 in arrays
    
    for i = 1:linelength
    
        x(:)=qu(t-1,i,:);
        y(:)=qu_y(t-1,i,:);
        first(:)=qu_first(t-1,i,:);
        second(:)=qu_second(t-1,i,:);
        
        qunew(t,i,:) = F*x;                     % operates on the qubit at each z          
        qunew_y(t,i,:) = F*y;
        qunew_first(t,i,:) = F*first;
        qunew_second(t,i,:) = F*second;
    end
    
    left = circshift(qunew,Kl);                 % determines which part moves left
    right = circshift(qunew,Kr);                % determines which part moves right
    left_y = circshift(qunew_y,Kl);
    right_y = circshift(qunew_y,Kr);
    left_first = circshift(qunew_first,Kl);
    right_first = circshift(qunew_first,Kr);
    left_second = circshift(qunew_second,Kl);
    right_second = circshift(qunew_second,Kr);
    
    qu(t,:,1) = left(t,:,1);                    % (a b) a part moves left
    qu(t,:,2) = right(t,:,2);                   % (a b) b part moves right
    qu_y(t,:,1) = left_y(t,:,1);
    qu_y(t,:,2) = right_y(t,:,2);
    qu_first(t,:,1) = left_first(t,:,1);
    qu_first(t,:,2) = right_first(t,:,2);
    qu_second(t,:,1) = left_second(t,:,1);
    qu_second(t,:,2) = right_second(t,:,2);
    
    prob(t,:) = (abs(qu(t,:,1)).^2 + abs(qu(t,:,2)).^2);
    prob_y(t,:) = (abs(qu_y(t,:,1)).^2 + abs(qu_y(t,:,2)).^2);
    prob_first(t,:) = (abs(qu_first(t,:,1)).^2 + abs(qu_first(t,:,2)).^2);
    prob_second(t,:) = (abs(qu_second(t,:,1)).^2 + abs(qu_second(t,:,2)).^2);
    
    if detector_everywhere == 1
        for ii = 1:linelength
            if t == detector_everywhere_iteration 
                if detect(prob(t,ii)) ==  1 
                    prob(t,:) = 0;
                    qu(t,:,1) = 0;
                    qu(t,:,2) = 0;
                    prob(t,ii) = 1;
                    qu(t,ii,1) = init_qubit(1,1);
                    qu(t,ii,2) = init_qubit(2,1);
                    break
                else
                    check = 0;     
                    prob(t,ii) = 0;
                    prob(t,:) = redistribute_prob(prob,t); 
                end
            end
        end
    end
        if detector_on == 1  
                for ii = 1:linelength
                    if detection_point(1,ii) == 1              % detect particle 
                        if prob(t,ii) ~= 0 && detect(prob(t,ii)) == 1 
                            qu = update_qu_detection(qu,t,ii,init_qubit);
                            prob(t,:) = 0;
                            prob(t,ii) = 1;
                        elseif prob(t,ii) ~= 0
                            j_sum = 0;
                            prob(t,ii) = 0;
                            sum_probs = 0;
                            for j = 1:linelength
                                %prob(t,j);
                                if prob(t,j) ~= 0 
                                    j_sum = j_sum + 1;
                                    prob;
                                    sum_probs = sum_probs + prob(t,j);
                                    update_qu_value(qu,t,linelength);
                                    %prob(t,:) = prob(t,:)*(1/sum_probs);
                                    qu = update_qu_nodetection(qu,t,ii,linelength);
                                end   
                            end
                            prob(t,:) = prob(t,:).*(1/sum_probs);
                            qu = update_qu_nodetection(qu,t,ii,linelength);
                            qu(t,ii,1) = 0;
                            qu(t,ii,2) = 0;
                        end
                            %prob(t,ii) = 0;
                    end 
                end
        elseif detect_superposition == 1
        	for ii = 1:linelength
            	if detection_point(1,ii) == 1              % detect particle 
                	if prob(t,ii) ~= 0 && detect(prob(t,ii)) == 1 
                    	qu = update_qu_detection(qu,t,ii,init_qubit);
                        	if choose_qubit == 1
                            	prob(t,:) = 0;
                            	prob(t,:) = prob(t,:) + prob_first(t,:);
                            	prob(t,ii) = prob(t,ii) + 1;
                            	qu(t,:) = (qu(t,:)+ update_qu_detection_first(qu_first,t,ii,init_qubit_first)+qu_second(t,:))/2;
                            else
                            	prob(t,:) = 0;
                             	prob(t,:) = prob(t,:) + prob_second(t,:);
                            	prob(t,ii) = prob(t,ii) + 1;
                            	qu(t,:) = (qu(t,:)+ update_qu_detection_second(qu_second,t,ii,init_qubit_second) + qu_second(t,:))/2;
                            end                                
                    elseif prob(t,ii) ~= 0
                    	j_sum = 0;
                    	prob(t,ii) = 0;
                    	sum_probs = 0;
                    	for j = 1:linelength
                        	%prob(t,j);
                        	if prob(t,j) ~= 0 
                            	j_sum = j_sum + 1;
                            	prob;
                            	sum_probs = sum_probs + prob(t,j);
                            	update_qu_value(qu,t,linelength);
                            	%prob(t,:) = prob(t,:)*(1/sum_probs);
                            	qu = update_qu_nodetection(qu,t,ii,linelength);
                            end   
                        end
                        prob(t,:) = prob(t,:).*(1/sum_probs);
                        qu = update_qu_nodetection(qu,t,ii,linelength);
                        qu(t,ii,1) = 0;
                        qu(t,ii,2) = 0;
                    end
                            %prob(t,ii) = 0;
                end 
            end                        
        end

    
    %line_num(qu,ntimes,linelength)
    prob(t,:) = (1/sum(prob(t,:)))*prob(t,:);    
    prob_sum = sum(prob,2);
     
    zprob(:,1) = prob(t,:).*z(t,:);
    z2(t,:) = z(t,:).*z(t,:);
    z2prob(:,1) = (prob(t,:).*z2(t,:));%.^2;
    
    mean_z(t,1) = sum(zprob,1);
    mean_z2(t,1) = sum(z2prob,1);
    
    stdev_z(t,1) = sqrt(abs(mean_z2(t,1) - mean_z(t,1)^2));
    
    vel_prob(:,1) = prob(t,:).*z(t,:) - prob(t-1,:).*z(t-1,:);  
    vel2_prob(:,1) = prob(t,:).*z2(t,:) - prob(t-1,:).*z2(t-1,:);
    
    vel_mn(t,1) = mean_z(t,1) - mean_z(t-1,1);
    vel2_mn(t,1) = mean_z2(t,1) - mean_z2(t-1,1);
    stdev_v(t,1) = sqrt(abs(vel2_mn(t,1) - vel_mn(t,1)^2));
    
    for xz = 1:linelength
        sum_x(t,1) = sum_x(t,1) + sym(prob(t,xz)*(xz-linecenter));
        sum_x2(t,1) = sum_x(t,1)^2;
    end
    
    uncertainty(t,1) = stdev_z(t,1).*stdev_v(t,1);
end
qu;


%for i = 2:ntimes;                       % method for finding average x 
  % for xz = 1:linelength;      
       %%%%end

  

fprintf(fid_meanZ,'        timestep           mean_z          mean_z2         st.dev-z         velocity       velocity^2         st.dev-v         prob chk        uncertainty\n');

for it = 1:ntimes
    
    wt(it,1) = 1;
    time(it,1) = it;
    qua = zeros(linelength,1);
    qub = zeros(linelength,1);
    
    qua(1:linelength,1) = qu(it,1:linelength,1);  %1,1:linelength,1
    qub(1:linelength,1) = qu(it,1:linelength,2);
    probLine(1:linelength,1) = prob(it,1:linelength);
    
    uncertainty_ = uncertainty(it,1);
    avg_x = sum_x(it,1);
    mean_z_t = mean_z(it,1);
    mean_z2_t = mean_z2(it,1);
    stdev_z_t = stdev_z(it,1);
    prob_sum_t = prob_sum(it,1);
    vel_t = vel_mn(it,1);
    vel2_t = vel2_mn(it,1);
    %stdev_v_t = std(vel_mn,wt);
    stdev_v_t = stdev_v(it,1);
    time_t = time(it,1);
    
    fprintf(fid_qubit,' %f ', real(qua));
        fprintf(fid_qubit,'     \n'); 
    fprintf(fid_qubit,' %f ', imag(qua));
        fprintf(fid_qubit,'     \n');
    fprintf(fid_qubit,' %f ', real(qub));
        fprintf(fid_qubit,'     \n'); 
    fprintf(fid_qubit,' %f ', imag(qub));
        fprintf(fid_qubit,'     \n');
        fprintf(fid_qubit,'     \n');
        
    fprintf(fid_probs,' %f ', probLine);
    fprintf(fid_probs,'     \n'); 
        
    fprintf(fid_meanZ,'%16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f \n', ...
        [time_t mean_z_t mean_z2_t stdev_z_t vel_t vel2_t stdev_v_t prob_sum_t uncertainty_]);
     
    counter_real_qua = nonzeros(real(qua)');
    counter_imag_qua = nonzeros(imag(qua)');
    counter_real_qub = nonzeros(real(qub)');
    counter_imag_qub = nonzeros(imag(qub)');
    real(qua);
    imag(qua);
    real(qub);
    imag(qub);
    
    
end


% figure(6);
% a = plot([1:length(counter_real_qua)],counter_real_qua);
% hold on
% b = plot([1:length(counter_imag_qub)],counter_imag_qub);
% hold off
% xlabel('z');
% ylabel('value');
% title('Real qua');
% 
% figure(8);
% a = plot([1:length(counter_imag_qua)],counter_imag_qua);
% xlabel('z');
% ylabel('value');
% title('imag qua');
% 
% figure(9);
% a = plot([1:length(counter_real_qub)],counter_real_qub);
% xlabel('z');
% ylabel('value');
% title('Real qub');
% 
% figure(5);
% a = plot([1:length(counter_imag_qub)],counter_imag_qub);
% hold on
% b = plot([1:length(counter_real_qub)],counter_real_qub);
% hold off
% xlabel('z');
% ylabel('value');
% title('imag qub');
% 
% figure(4);
% a = plot([1:length(real(qua))],real(qua));
% hold on
% b = plot([1:length(imag(qua))],(imag(qua)));
% c = plot([1:length(real(qub))],(real(qub)));
% d = plot([1:length(imag(qub))],(imag(qub))); 
% legend([a,b,c,d], 'real(a)','imag(a)','real(b)','imag(b)');
% hold off
% xlabel('z');
% ylabel('value');
% title('zeros qua qub comparison');
% 
% figure(4);
% a = plot([1:length(counter_real_qua)],counter_real_qua);
% hold on
% b = plot([1:length(counter_imag_qua)],counter_imag_qua);
% c = plot([1:length(counter_real_qub)],counter_real_qub);
% d = plot([1:length(counter_imag_qub)],counter_imag_qub); 
% legend([a,b,c,d], 'real(a)','imag(a)','real(b)','imag(b)');
% hold off
% xlabel('z');
% ylabel('value of number');
% title('No zeros qua qub comparison');
% 
% figure(6);
% plot(time,vel_mn);
% title('vel vs. time');
% xlabel('t');
% ylabel('vel');
% axis([0 ntimes -1.2 1.2]);
% 
Image_zt(prob,1,'parula');   
% 
% 
% Image_zt(prob,1,'parula');
% Image_zt(prob_first,2,'parula');
% Image_zt(prob_second,3,'parula');
% 
% Image_zt(prob_y,2,'parula')
% 
% figure(5)
% plot(time,uncertainty);
% title('Uncertainty');
% xlabel('Iteration');
% ylabel('uncertainty');

%--------------------------------------------------------------
% VISUALIZATION
% Copyright FauxPlay 2019
function dummy = Image_zt(Universe,figno,cmap)

figure(figno)
colormap(cmap);
imagesc(Universe,[0 1]);
colorbar;
title('probability density');
set(gca,'FontSize',10)
drawnow

% frame = getframe(gcf);
% writeVideo(v,frame);

end

function output = detect(qubit)
    detection_rate  = rand;
    if detection_rate < qubit 
        output = 1; 
    else 
        output = 0; 
    end
end 


function output = choose_qubit
    choose = rand;
    if choose < 0.5
        output = 1;
    else 
        output = 0;
    end 
end

function output = update_qu_detection(qu,t,ii,init_qubit)
    qu(t,:,1) = 0;
    qu(t,ii,1) = init_qubit(1,1);
    qu(t,:,2) = 0;
    qu(t,ii,2) = init_qubit(2,1);
    output = qu;
end

function output = update_qu_detection_first(qu_first,t,ii,init_qubit_first)
    qu_first(t,:,1) = 0;
    qu_first(t,ii,1) = init_qubit_first(1,1);
    qu_first(t,:,2) = 0;
    qu_first(t,ii,2) = init_qubit_first(2,1);
    output = qu_first(t,:);
end

function output = update_qu_detection_second(qu_second,t,ii,init_qubit_second)
    qu_second(t,:,1) = 0;
    qu_second(t,ii,1) = init_qubit_second(1,1);
    qu_second(t,:,2) = 0;
    qu_second(t,ii,2) = init_qubit_second(2,1);
    output = qu_second(t,:);
end

function output = line_num(qu,ntimes,linelength)
sum_nonzero = 0;
for f = 2:ntimes
    for g = 1:linelength
        if qu(f,g,1) ~= 0
            sum_nonzero = sum_nonzero + 1;
        end
    end
end
output = sum_nonzero;
end

function output = update_qu_nodetection(qu,t,ii,linelength)
qu(t,:,1);
qu1_update(t,:,1) = qu(t,:,1);
qu1_update(t,ii,1) = 0;
qu2_update(t,:,2) = qu(t,:,2);
qu2_update(t,ii,2) = 0;
factor_qu = sum(abs(qu1_update(t,:,1))+abs(qu2_update(t,:,2)));
divide = update_qu_value(qu,t,linelength);
for k = 1:linelength
    if qu(t,k,1) ~= 0
        if qu(t,k,1) > 0 
            qu(t,k,1) = qu(t,k,1);%*(2*normalize/factor_qu);% + qu(t,ii,1)/divide;
%         else
%             qu(t,k,1) = qu(t,k,1) - qu(t,ii,1)/divide;
        end 
    end
    if qu(t,k,2) ~= 0
        if qu(t,k,1) > 0 
            qu(t,k,2) = qu(t,k,2);%*(2*normalize/factor_qu);% + qu(t,ii,2)/divide;
%         else
%             qu(t,k,2) = qu(t,k,2) - qu(t,ii,2)/divide;
        end 
    end
end
output = qu;
end

function output = update_qu_value(qu,t,linelength)
output = 0;
for i = 1:linelength
    test = qu(t,i,1);
    if qu(t,i,1) ~= 0
       output = output + 1;
    else 
       output = output;
    end
end
end

function output = redistribute_prob(prob,t)
x = 1/sum(prob(t,:));
prob(t,:) = x*prob(t,:);
output = prob(t,:);
end
