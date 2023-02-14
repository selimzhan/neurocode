clear all
close all
clc
addpath(genpath('Z:\Lab\Lesion_Dannce\L2_2\'))

%% load labels from the same session
session1 = load('Z:\Lab\Lesion_Dannce\L2_2\20210331\20220713_122437_Label3D_dannce.mat');
%session2 = load('Z:\Lab\Lesion_Dannce\L2_2\20210331\20220712_141243_Label3D_dannce.mat');
session3 = load('Z:\Lab\Lesion_Dannce\L2_2\20210331\DANNCE\predict_new_2\save_data_AVG.mat');
%session4 = load('C:\Users\selim\Desktop\2021_8_4_2000_BehaviorToneWater\dannce\predict12\save_data_AVG');

%% Get the skeleton 
%skeleton = load('C:\Users\selim\Desktop\skeletons\rat14.mat');
skeleton = load('Z:\Lab\Lesion_Dannce\Label4D\Label4D\skeletons\jesse_skeleton.mat');
npts = 20;

%% extract desired frames
dframes = session1.labelData{1,1}.data_frame;

%% extract desired coordinatess
pts1 = session1.labelData{1,1}.data_3d;
%pts2 = session2.labelData{1,1}.data_3d;
pts3 = session3.pred(dframes,:,:);
%pts4 = session4.pred(dframes,:,:);

%% initialize data
nframes = size(dframes, 2);
coordinates1 = zeros(npts,3,nframes);
%coordinates2 = zeros(npts,3,nframes);
dist = zeros(nframes,14,3);

%% calculate the euclidean distance for each point from the same frame and plot labels 
allpts = npts*3;

for n=1:nframes
    coordinates1(:,1,n) = pts1(n,1:3:allpts)';
    coordinates1(:,2,n) = pts1(n,2:3:allpts)'; 
    coordinates1(:,3,n) = pts1(n,3:3:allpts)';
    %coordinates2(:,1,n) = pts2(n,1:3:allpts)';
    %coordinates2(:,2,n) = pts2(n,2:3:allpts)';
    %coordinates2(:,3,n) = pts2(n,3:3:allpts)';
    
    figure(n)
    plot3(squeeze(coordinates1(:,1,n)), squeeze(coordinates1(:,2,n)), squeeze(coordinates1(:,3,n)), 'or');
    %plot3(squeeze(coordinates2(:,1,n)), squeeze(coordinates2(:,2,n)), squeeze(coordinates2(:,3,n)), 'og');
    plot3(squeeze(pts3(n,1,:)), squeeze(pts3(n,2,:)), squeeze(pts3(n,3,:)), 'ob');
    %plot3(squeeze(pts4(n,1,:)), squeeze(pts4(n,2,:)), squeeze(pts4(n,3,:)), 'ob');
    for p=1:npts
        %dist(n,p,1) = norm(coordinates1(p,:,n)-coordinates2(p,:,n));
        dist(n,p,2) = norm(coordinates1(p,:,n)-pts3(n,:,p));
        %dist(n,p,3) = norm(coordinates2(p,:,n)-pts3(n,:,p));
        idx1 = skeleton.joints_idx(p,1);
        idx2 = skeleton.joints_idx(p,2);
        plot3([coordinates1(idx1,1,n) coordinates1(idx2,1,n)], [coordinates1(idx1,2,n) coordinates1(idx2,2,n)], [coordinates1(idx1,3,n) coordinates1(idx2,3,n)],'r')
        hold on
        %plot3([coordinates2(idx1,1,n) coordinates2(idx2,1,n)], [coordinates2(idx1,2,n) coordinates2(idx2,2,n)], [coordinates2(idx1,3,n) coordinates2(idx2,3,n)], 'g')
        plot3([pts3(n,1,idx1) pts3(n,1,idx2)], [pts3(n,2,idx1) pts3(n,2,idx2)], [pts3(n,3,idx1) pts3(n,3,idx2)], 'b')
        %plot3([pts4(n,1,idx1) pts4(n,1,idx2)], [pts4(n,2,idx1) pts4(n,2,idx2)], [pts4(n,3,idx1) pts4(n,3,idx2)], 'b')
        %h1 = plot3([coordinates1(idx1,1,n) coordinates1(idx2,1,n)], [coordinates1(idx1,2,n) coordinates1(idx2,2,n)], [coordinates1(idx1,3,n) coordinates1(idx2,3,n)]);
        %h2 = plot3([pts3(n,1,idx1) pts3(n,1,idx2)], [pts3(n,2,idx1) pts3(n,2,idx2)], [pts3(n,3,idx1) pts3(n,3,idx2)], 'g');
        %h3 = plot3([pts4(n,1,idx1) pts4(n,1,idx2)], [pts4(n,2,idx1) pts4(n,2,idx2)], [pts4(n,3,idx1) pts4(n,3,idx2)], 'b');
        %colors = {[237/255 85/255 59/255], [60/255 174/255 163/255], [32/255 99/255	155/255]};
        %[h1.Color, h2.Color, h3.Color] = colors{:};
        set(gcf, 'renderer', 'painters')
    end
    hold off
end 

dist_k_s_1_pr_11 = dist(:,:,1);
%dist_d_s_1_pr_11 = dist(:,:,2);
dist_d_k_1_pr_11 = dist(:,:,3);

%%
save('distances_between_users_L1_4.mat','dist_k_s_1_pr_11','dist_d_s_1_pr_11','dist_d_k_1_pr_11')
%% sTaTIstIcaL meAsuREs (not very rigorous)
mean_arr = mean(dist);
mean_2d = cat(1, mean_arr(1,:,1), mean_arr(1,:,2), mean_arr(1,:,3));
mean_dannce_user = (mean_2d(2,:) + mean_2d(3,:))/2;
new_mean_2d = cat(1, mean_2d(1,:), mean_dannce_user);
mean_over_points = mean(new_mean_2d,2);
%points_std = num2cell(mean_arr);
%body_parts_std = [skeleton.joint_names'; points_std(:,:,1)]
%body_parts_std_next = [body_parts_std; points_std(:,:,2)]
%body_parts_std_final = [body_parts_std_next; points_std(:,:,3)]

%% difference bars across joints
bar(new_mean_2d,2)