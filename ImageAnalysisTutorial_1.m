% Image analysis exercise: 
% 1) play with importing images in matlab
% 2) play with measuring intensities 
% play with segmenting out bright objects from an image

%% Select directory on your computer containing images

folder_name = uigetdir; % opens a GUI to let you navigate to your chosen folder of images
% QUESTION 1: In your workspace, check out the value of the variable
% "folder_name". What did uigetdir ("user interface get directory" do?


cd(folder_name); % makes that folder matlab's current working folder
% QUESTION 2: Check your "current folder" window within matlab. Is it in the right
% place?  Do you see your list of images?


%% Now, we'll import an image into a matlab array

% In the line below, replace 'name_of_your_image.tif' with your chosen image?

myImage  = imread('colony_growth_25.tif'); % import your image
% QUESTION 3: In your workspace, open and check out "myImage".  Think about
% what an image is. If you move to the image of your Matlab array, can you
% find the bright aster from your image?

%% Display image visually
figure
imshow(myImage, []) % Display the image. Adding the [] scales how the image displays, 
% so that the max value is white and the min is black, and everything in between is gray
% QUESTION 4: Try removing the '[]'. Do things look different?

%% Measuring intensity

% Measure the average intensity in the top left corner of the image.
% Matlab arrays define (0,0) in the top left corner of an array. Since
% myImage is an array of pixel values, how would you find the average
% intensity in a 10 pixel x 10 pixel square at the top left?

sumTopLeftCorner = sum(sum(myImage(1:10,1:10) )); % sum up intensity values across 10 rows and 10 columns
averageTopLeftCorner = sumTopLeftCorner / 100; % divide by 100 boxes to get average intensity

% QUESTION 5: What's the average intensity in a box at the center of the
% image that's 30 boxes tall and 20 boxes wide?

%% Playing with segmentation

% To segment out bright things in an image, we'll set a threshold - and
% throw out all the pixels lower than that threshold by setting them to zero

% First, let's get some clues about what our threshold should be
figure
imhist(myImage, 500); %histogram of image data with 500 bins

% Based on our histogram, can you guess a reasonable threshold value to
% try?

%% Threshold

% Side note: what is thresholding, exactly?
% QUESTION 6: copy and paste the two lines below into your command window.
% Open up the array 'a' and 'aThresh' and compare. 'aThresh' is a "binary
% mask of a."  Does that make sense?
a = [2,4,1,6,4,3; 5,2,7,1,2,4]; % make an example array
aThresh = a > 3; % aThresh will have a 1 in all boxes where a>3, and a 0 where a<0

% Now let's do our image
threshold = 100; % set threshold to your first guess (delete '100')
% This is up to you - maybe you want just the aster core, or most of the
% aster, or whatever

imThresh = myImage > threshold; % this does the thresholding

figure
imshow(imThresh); % check our segmentation. 
% If it's not what you wanted, just change the value of threshold and repeat this section.

%% Another possibly useful thing: measure the area of the segmented region

% Measure the area of the segmented region
area = sum(sum(imThresh)); % because our image is binary (0's and 1's), this sum will tell us the AREA
% of the segmented region

% QUESTION 7: Does this area roughly make sense to you?  What fraction is
% it of the total image?

%% Another possibly useful thing

% Extract the concentration profile of whatever's in your image

% Open up your original image in Matlab's imtool
imtool(myImage)

% QUESTION 8: what pixel value represents the center of your aster?

%% QUESTION 9: given that center, can you pull out a vector of intensity
% values that starts at the center of your aster and goes out to one edge?
% In other words, this might look like the "exponential distributions" we
% saw in Rachel's data. Something like:

% PS it's possible you have to close the imtool?

initialProfile = myImage(200, 200:400); % this might work if your aster center were (300, 300)
% and (200, 400) was near the right edge of your image.  Does that make
% sense?

% Check out your initialProfile in your matlab workspace to see if it looks
% right, or plot it
xVector = (1:length(initialProfile));
figure
plot(xVector, initialProfile)
xlabel('position (pixels)')
ylabel('intensity (AU)')









