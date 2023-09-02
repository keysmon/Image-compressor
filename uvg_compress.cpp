/* uvg_compress.cpp

   Starter code for Assignment 3 (in C++). This program
    - Reads an input image in BMP format
     (Using bitmap_image.hpp, originally from 
      http://partow.net/programming/bitmap/index.html)
    - Transforms the image from RGB to YCbCr (i.e. "YUV").
    - Downscales the Cb and Cr planes by a factor of two
      (producing the same resolution that would result
       from 4:2:0 subsampling, but using interpolation
       instead of ignoring some samples)
    - Writes each colour plane (Y, then Cb, then Cr)
      in 8 bits per sample to the output file.

   B. Bird - 2023-07-03
*/
#include <cmath>
#include <math.h> 
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include <cstdint>
#include "output_stream.hpp"
#include "bitmap_image.hpp"
#include "uvg_common.hpp"
using namespace std;
double block_size = 8;
#define pi 3.14

unsigned int height;
unsigned int width;

vector<vector<double>> quantization_matrix_lumin = {
{16,11,10,16,24,40,51,61},
{12,12,14,19,26,58,60,55},
{14,13,16,24,40,57,69,56},
{14,17,22,29,51,87,80,62},
{18,22,37,56,68,109,103,77},
{24,35,55,64,81,104,113,92},
{49,64,78,87,103,121,120,101},
{72,92,95,98,112,100,103,99}};

vector<vector<double>> quantization_matrix_chrom = {
{17,18,24,47,99,99,99,99},
{18,21,26,66,99,99,99,99},
{24,26,56,99,99,99,99,99},
{47,66,99,99,99,99,99,99},
{99,99,99,99,99,99,99,99},
{99,99,99,99,99,99,99,99},
{99,99,99,99,99,99,99,99},
{99,99,99,99,99,99,99,99}};


vector<vector<int>> zig_zag_route = {
{0,0},
{0,1},{1,0},
{0,2},{1,1},{2,0},
{0,3},{2,1},{1,2},{3,0},
{0,4},{3,1},{2,2},{1,3},{4,0},
{0,5},{4,1},{3,2},{2,3},{1,4},{5,0},
{0,6},{5,1},{4,2},{3,3},{4,2},{1,5},{6,0},
{0,7},{6,1},{2,5},{3,4},{4,3},{5,2},{1,6},{7,0},
{1,7},{2,6},{3,5},{4,4},{5,3},{6,2},{7,1},
{2,7},{3,6},{4,5},{5,4},{6,3},{7,2},
{3,7},{4,6},{5,5},{6,4},{7,3},
{4,7},{5,6},{6,5},{7,4},
{5,7},{6,6},{7,5},
{6,7},{7,6},
{7,7}};


vector<vector<int>> zig_zag_coords(int height, int width,int frequent_size){
    vector<vector<int>> coords = {};
    int num_block_height = int(ceil(double(height)/block_size));
    int num_block_width = int(ceil(double(width)/block_size));
    
    for (int y = 0;y < num_block_height ; y++){
        for (int x = 0 ; x < num_block_width ; x++){
            for (int i = 0 ; i< frequent_size;i++){
                int y_cord = zig_zag_route[i][0] + y * block_size;
                int x_cord = zig_zag_route[i][1] + x * block_size;
                if (y_cord < height && x_cord < width){
                    vector<int> temp = {y_cord,x_cord};
                    coords.push_back(temp);
                }
            }
        }
    }
    return coords;

}


// Function to find discrete cosine transform and print itv
// implementaion inpsired by https://www.geeksforgeeks.org/discrete-cosine-transform-algorithm-program/
vector<vector<double>> dctTransform(vector<vector<double>> matrix)
{   
    #define pi 3.14
    int i, j, k, l;
    int len = 8;
    // dct will store the discrete cosine transform
    vector<vector<double>> dct = create_2d_vector<double>(block_size,block_size);
    double ci, cj, dct1, sum;
    for (i = 0; i < len; i++) {
        for (j = 0; j < len; j++) {
            
            // ci and cj depends on frequency as well as
            // number of row and columns of specified matrix
            if (i == 0)
                ci = 1 / sqrt(block_size);
            else
                ci = sqrt(2) / sqrt(block_size);
            if (j == 0)
                cj = 1 / sqrt(block_size);
            else
                cj = sqrt(2) / sqrt(block_size);

            // sum will temporarily store the sum of
            // cosine signals
            sum = 0;
            for (k = 0; k < block_size; k++) {
                for (l = 0; l < block_size; l++) {
                    dct1 = matrix[k][l] *
                           cos((2 * k + 1) * i * pi / (2 * block_size)) *
                           cos((2 * l + 1) * j * pi / (2 * block_size));
                    sum = sum + dct1;
                }
            }
            dct[i][j] = ci * cj * sum;
        }
    }
    return dct;
}



vector<vector<unsigned char>> DCT (vector<vector<unsigned char>> image,vector<vector<double>> quanti_matrix){
    int num_block_height = int(ceil(double(image.size())/block_size));
    int num_block_width = int(ceil(double(image[0].size())/block_size));
    int padded_height = num_block_height * block_size;
    int padded_width = num_block_width * block_size;

    // add padding to the image so that the size is multiples of block_size
    vector<vector<double>> padded_image(padded_height, vector<double>(padded_width, 0));
    //vector<vector<double>> padded_image(padded_height, vector<double>(padded_width, 0));
    for(int y = 0; y < image.size(); y++){
        for (int x = 0; x < image[0].size(); x++){
            padded_image[y][x] =  image[y][x] - 128;
        }
    }
    
    for (int y = 0; y < image.size();y++){
        for (int x = image[0].size(); x < padded_image[0].size();x++){
            padded_image[y][x] = image[y][image[0].size()-1] - 128;
        }
    }
    
    
    for (int x = 0; x < padded_image[0].size();x++){
        for (int y = image.size(); y < padded_image.size();y++){
            padded_image[y][x] = image[image.size()-1][x] - 128;
            
        }
    }
    auto block = create_2d_vector<double>(int(block_size),int(block_size));

    for (int col_idx = 0; col_idx < num_block_height ; col_idx++){
        int col_idx_start = col_idx * block_size;
        
        for (int row_idx = 0; row_idx < num_block_width  ; row_idx++){
            int row_idx_start = row_idx * block_size;
            
            // extract the contents in the block
            for (int i = 0;i<block_size;i++){
                for (int j =0;j<block_size;j++){
                    block[i][j] = padded_image[col_idx_start + i][row_idx_start + j];
                }
            }
            /*
            for (auto i : block){
                for (auto j : i){
                    cout << int(j) << " ";
                }cout << endl;
            }cout << endl;*/
            // DCT transform
            auto dct_block = dctTransform(block);
            
            // Quantization transform
            for (size_t i = 0; i < block_size; i++) {
                for (size_t j = 0; j < block_size; j++) {
                    dct_block[i][j] = round(dct_block[i][j] / quanti_matrix[i][j]);
                }
            }
            // copy blocks back to image and add 128 offset to remove negative sign
            for (int i = 0;i<block_size;i++){
                for (int j=0;j<block_size;j++){
                    padded_image[col_idx_start+i][row_idx_start+j] = dct_block[i][j]+128;
                }
            }
        }
    }
    for (int i = 0; i<image.size();i++){
        for (int j =0;j<image[0].size();j++){
            image[i][j] = padded_image[i][j];
        }
    }
   

    return image;
    
}

//A simple downscaling algorithm using averaging.
std::vector<std::vector<unsigned char> > scale_down(std::vector<std::vector<unsigned char> > source_image, unsigned int source_width, unsigned int source_height, int factor){

    unsigned int scaled_height = (source_height+factor-1)/factor;
    unsigned int scaled_width = (source_width+factor-1)/factor;

    //Note that create_2d_vector automatically initializes the array to all-zero
    auto sums = create_2d_vector<unsigned int>(scaled_height,scaled_width);
    auto counts = create_2d_vector<unsigned int>(scaled_height,scaled_width);

    for(unsigned int y = 0; y < source_height; y++)
        for (unsigned int x = 0; x < source_width; x++){
            sums.at(y/factor).at(x/factor) += source_image.at(y).at(x);
            counts.at(y/factor).at(x/factor)++;
        }
    auto result = create_2d_vector<unsigned char>(scaled_height,scaled_width);
    for(unsigned int y = 0; y < scaled_height; y++)
        for (unsigned int x = 0; x < scaled_width; x++)
            result.at(y).at(x) = (unsigned char)((sums.at(y).at(x)+0.5)/counts.at(y).at(x));
    return result;
}



int main(int argc, char** argv){
    


    if (argc < 4){
        std::cerr << "Usage: " << argv[0] << " <low/medium/high> <input BMP> <output file>" << std::endl;
        return 1;
    }
    std::string quality{argv[1]};
    std::string input_filename {argv[2]};
    std::string output_filename {argv[3]};
    
    bitmap_image input_image {input_filename};
    height = input_image.height();
    width = input_image.width();
    //Read the entire image into a 2d array of PixelRGB objects 
    //(Notice that height is the outer dimension, so the pixel at coordinates (x,y) 
    // must be accessed as imageRGB.at(y).at(x)).
    std::vector<std::vector<PixelYCbCr>> imageYCbCr = create_2d_vector<PixelYCbCr>(height,width);
    
    for(unsigned int y = 0; y < height; y++){
        for (unsigned int x = 0; x < width; x++){
            auto [r,g,b] = input_image.get_pixel(x,y);
            PixelRGB rgb_pixel {r,g,b};
            imageYCbCr.at(y).at(x) = rgb_pixel.to_ycbcr();            
        }
    }
    

    
    std::ofstream output_file{output_filename,std::ios::binary};
    OutputBitStream output_stream {output_file};

    //Placeholder: Use a simple bitstream containing the height/width (in 32 bits each)
    //followed by the entire set of values in each colour plane (in row major order).
    int frequent_size = 0;

    //cout << "Height - " << height << "  Width - " << width<<endl;
    output_stream.push_u32(height);
    output_stream.push_u32(width);
    
    if (quality == "low")
        frequent_size = 1;
    else if (quality == "medium")
        frequent_size = 6;
    else if (quality == "high")
        frequent_size = 32;
    else if (quality == "full")
        frequent_size = 64;
    output_stream.push_byte(frequent_size);

    
    //Extract the Y plane into its own array 
    auto Y = create_2d_vector<unsigned char>(height,width);
    for(unsigned int y = 0; y < height; y++)
        for (unsigned int x = 0; x < width; x++)
            Y.at(y).at(x) = imageYCbCr.at(y).at(x).Y;

    //Extract the Cb plane into its own array 
    auto Cb = create_2d_vector<unsigned char>(height,width);
    for(unsigned int y = 0; y < height; y++)
        for (unsigned int x = 0; x < width; x++)
            Cb.at(y).at(x) = imageYCbCr.at(y).at(x).Cb;
    auto Cb_scaled = scale_down(Cb,width,height,2);

    //Extract the Cr plane into its own array 
    
    auto Cr = create_2d_vector<unsigned char>(height,width);
    for(unsigned int y = 0; y < height; y++)
        for (unsigned int x = 0; x < width; x++)
            Cr.at(y).at(x) = imageYCbCr.at(y).at(x).Cr;
    auto Cr_scaled = scale_down(Cr,width,height,2);

    auto DCT_Y = DCT(Y,quantization_matrix_lumin);
    auto DCT_Cb = DCT(Cb_scaled,quantization_matrix_chrom);
    auto DCT_Cr = DCT(Cr_scaled,quantization_matrix_chrom);
    

    /*
    //Write the Y values 
    for(unsigned int y = 0; y < height; y++)
        for (unsigned int x = 0; x < width; x++)
            if ( (x%8 < frequent_size) && (y%8 < frequent_size) ){
                output_stream.push_byte(DCT_Y.at(y).at(x));
            }
    //Write the Cb values 
    for(unsigned int y = 0; y < (height+1)/2; y++)
        for (unsigned int x = 0; x < (width+1)/2; x++)
            if ( (x%8 < frequent_size) && (y%8 < frequent_size) ){

                output_stream.push_byte(DCT_Cb.at(y).at(x));
            }
    //Write the Cr values 
    for(unsigned int y = 0; y < (height+1)/2; y++)
        for (unsigned int x = 0; x < (width+1)/2; x++)
            if ( (x%8 < frequent_size) && (y%8 < frequent_size) ){
                output_stream.push_byte(DCT_Cr.at(y).at(x));
            }
    
    */
    auto Y_cords = zig_zag_coords(height,width,frequent_size);
    auto Cb_cords = zig_zag_coords((height+1)/2,(width+1)/2,frequent_size);
    auto Cr_cords = zig_zag_coords((height+1)/2,(width+1)/2,frequent_size);

    for (int i = 0;i<Y_cords.size();i++){
        output_stream.push_byte(DCT_Y.at(Y_cords[i][0]).at(Y_cords[i][1]));
    }
    for (int i = 0;i<Cb_cords.size();i++){
        output_stream.push_byte(DCT_Cb.at(Cb_cords[i][0]).at(Cb_cords[i][1]));
    }

    for (int i = 0;i<Cr_cords.size();i++){
        output_stream.push_byte(DCT_Cr.at(Cr_cords[i][0]).at(Cr_cords[i][1]));
    }

    output_stream.flush_to_byte();
    output_file.close();

    return 0;
}