/* uvg_decompress.cpp

   Starter code for Assignment 3 (in C++). This program
    - Reads a height/width value from the input file
    - Reads YCbCr data from the file, with the Y plane
      in full w x h resolution and the other two planes
      in half resolution.
    - Upscales the Cb and Cr planes to full resolution and
      transforms them to RGB.
    - Writes the result as a BMP image
     (Using bitmap_image.hpp, originally from 
      http://partow.net/programming/bitmap/index.html)

   B. Bird - 2023-07-03
*/

#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <cassert>
#include <cstdint>
#include "input_stream.hpp"
#include "bitmap_image.hpp"
#include "uvg_common.hpp"
using namespace std;
double block_size = 8;



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
/*
vector<vector<double>> quantization_matrix_lumin = {
{10,10,10,10,10,10,10,10},
{10,10,10,10,10,10,10,10},
{10,10,10,10,10,10,10,10},
{10,10,10,10,10,10,10,10},
{10,10,10,10,10,10,10,10},
{10,10,10,10,10,10,10,10},
{10,10,10,10,10,10,10,10},
{10,10,10,10,10,10,10,10}};

vector<vector<double>> quantization_matrix_chrom = {
{10,10,10,10,10,10,10,10},
{10,10,10,10,10,10,10,10},
{10,10,10,10,10,10,10,10},
{10,10,10,10,10,10,10,10},
{10,10,10,10,10,10,10,10},
{10,10,10,10,10,10,10,10},
{10,10,10,10,10,10,10,10},
{10,10,10,10,10,10,10,10}};
*/
// Function to find discrete cosine transform and print it



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

vector<vector<double>> inverse_dctTransform(vector<vector<double>> matrix)
{   
    #define pi 3.14
    int i, j, x, y;
   
    // dct will store the discrete cosine transform
    vector<vector<double>> dct = create_2d_vector<double>(block_size,block_size);
    double ci, cj, dct1, sum;
    for (x = 0;  x < block_size; x++) {
        for (y = 0; y < block_size; y++) {
            // ci and cj depends on frequency as well as
            // number of row and columns of specified matrix
            
            // sum will temporarily store the sum of
            // cosine signals
            sum = 0;
            for (i = 0; i < block_size; i++) {
                for (j = 0; j < block_size; j++) {
                    if (i == 0)
                        ci = 1 / sqrt(2);
                    else
                        ci = 1;
                    if (j == 0)
                        cj = 1 / sqrt(2);
                    else
                        cj = 1;
                    dct1 = matrix[i][j] *
                           cos((2 * x + 1) * i * pi / (2 * block_size)) *
                           cos((2 * y + 1) * j * pi / (2 * block_size)) /block_size;
                    sum = sum + dct1;
                }
            }
            dct[x][y] = ci * cj * sum;
        }
    }
    /*
    for (i = 0; i < block_size; i++) {
        for (j = 0; j < block_size; j++) {
            printf("%f\t", dct[i][j]);
        }
        printf("\n");
    }
    cout << "finsihed " <<endl; 
    */
    return dct;
}

vector<vector<unsigned char>> inverse_DCT (vector<vector<unsigned char>> image,vector<vector<double>> quanti_matrix){

    int num_block_height = int(ceil(double(image.size())/block_size));
    int num_block_width = int(ceil(double(image[0].size())/block_size));
    int padded_height = num_block_height * block_size;
    int padded_width = num_block_width * block_size;
    //cout << padded_height <<endl;
    //cout << padded_width <<endl;
    
   
    
    // add padding to the image so that the size is multiples of block_size
    vector<vector<double>> padded_image(padded_height, vector<double>(padded_width, 0));
    for(int y = 0; y < image.size(); y++){
        for (int x = 0; x < image[0].size(); x++){
            padded_image[y][x] =  image[y][x] -128;
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
                    block[i][j] = padded_image[col_idx_start + i][row_idx_start+j];
                }
            }

            /*
            for (int i = 0;i < padded_image.size();i++){
                for (int j = 0;j<padded_image[0].size();j++){
                    cout << padded_image[i][j] << " " ;
                }cout << endl<<endl;
            }*/

            // Quantization transform
            for (size_t i = 0; i < block_size; i++) {
                for (size_t j = 0; j < block_size; j++) {
                    block[i][j] = block[i][j] * quanti_matrix[i][j];
                }
            }
            
            /*
            for (int i = 0;i < block.size();i++){
                for (int j = 0;j<block[0].size();j++){
                    cout << block[i][j] << " " ;
                }cout << endl;
            }cout << endl;
            */

            // DCT transform
            auto dct_block = inverse_dctTransform(block);
            /*
            for (auto i : dct_block){
                for (auto j : i){
                    cout << int(j) << " ";
                }cout << endl;
            }cout << endl;*/

            /*
            for (int i = 0;i < padded_image.size();i++){
                for (int j = 0;j<padded_image[0].size();j++){
                    cout << padded_image[i][j] << " " ;
                }cout << endl<<endl;
            }*/

            for (int i = 0;i<block_size;i++){
                for (int j=0;j<block_size;j++){
                    padded_image[col_idx_start+i][row_idx_start+j] = dct_block[i][j];
                }
            }
        }
    }
    /*
    for (int i = 0;i < padded_image.size();i++){
        for (int j = 0;j<padded_image[0].size();j++){
            cout << padded_image[i][j] << " " ;
        }cout << endl<<endl;
    }*/
    
    for (int i = 0; i<image.size();i++){
        for (int j =0;j<image[0].size();j++){
            image[i][j] = round_and_clamp_to_char(padded_image[i][j] + 128);
        }
    }

    return image;
}




int main(int argc, char** argv){
   
    
    if (argc < 3){
        std::cerr << "Usage: " << argv[0] << " <input file> <output BMP>" << std::endl;
        return 1;
    }
    std::string input_filename {argv[1]};
    std::string output_filename {argv[2]};

    std::ifstream input_file{input_filename,std::ios::binary};
    InputBitStream input_stream {input_file};

    height = input_stream.read_u32();
    width = input_stream.read_u32();
    int frequent_size = input_stream.read_byte();
    //cout << "Frequent_size - " << frequent_size <<endl;
    //cout << "Height - " << height << "  Width - " << width<<endl;

    /*
    auto Y = create_2d_vector<unsigned char>(height,width);
    auto Cb_scaled = create_2d_vector<unsigned char>((height+1)/2,(width+1)/2);
    auto Cr_scaled = create_2d_vector<unsigned char>((height+1)/2,(width+1)/2);
    */

    vector<vector<unsigned char>> Y(height, vector<unsigned char>(width, 128));
    vector<vector<unsigned char>> Cb_scaled((height+1)/2, vector<unsigned char>((width+1)/2, 128));
    vector<vector<unsigned char>> Cr_scaled((height+1)/2, vector<unsigned char>((width+1)/2, 128));
    
    auto Y_cords = zig_zag_coords(height,width,frequent_size);
    auto Cb_cords = zig_zag_coords((height+1)/2,(width+1)/2,frequent_size);
    auto Cr_cords = zig_zag_coords((height+1)/2,(width+1)/2,frequent_size);
    int x1, x2;
    for (int i=0; i < Y_cords.size();i++){
        x1 = Y_cords[i][0];
        x2 = Y_cords[i][1];
        Y.at(x1).at(x2) = input_stream.read_byte();
    }
    for (int i=0; i < Cb_cords.size();i++){
        x1 = Cb_cords[i][0];
        x2 = Cb_cords[i][1];
        Cb_scaled.at(x1).at(x2) = input_stream.read_byte();
    }
    for (int i=0; i < Cr_cords.size();i++){
        x1 = Cr_cords[i][0];
        x2 = Cr_cords[i][1];
        Cr_scaled.at(x1).at(x2) = input_stream.read_byte();
    }


    /*
    for (unsigned int y = 0; y < height; y++){
        for (unsigned int x = 0; x < width; x++){
            if ( (x%8 < frequent_size) && (y%8 < frequent_size) ){
                Y.at(y).at(x) = input_stream.read_byte();
            }else{
                Y.at(y).at(x) = 128;
            }
        }
    }
    for (unsigned int y = 0; y < (height+1)/2; y++){
         for (unsigned int x = 0; x < (width+1)/2; x++){
             if ( (x%8 < frequent_size) && (y%8 < frequent_size) ){
                Cb_scaled.at(y).at(x) = input_stream.read_byte();
            }else{
                Cb_scaled.at(y).at(x) = 128;
            }
        }
    }
    for (unsigned int y = 0; y < (height+1)/2; y++){
        for (unsigned int x = 0; x < (width+1)/2; x++){
            if ( (x%8 < frequent_size) && (y%8 < frequent_size) ){
                Cr_scaled.at(y).at(x) = input_stream.read_byte();
            }else{
                Cr_scaled.at(y).at(x) = 128;
            }
        }
    }
    */

    auto Y_inverse = inverse_DCT(Y,quantization_matrix_lumin);
    auto Cb_scaled_inverse = inverse_DCT(Cb_scaled,quantization_matrix_chrom);
    auto Cr_scaled_inverse = inverse_DCT(Cr_scaled,quantization_matrix_chrom);
    

    auto imageYCbCr = create_2d_vector<PixelYCbCr>(height,width);
    for (unsigned int y = 0; y < height; y++){
        for (unsigned int x = 0; x < width; x++){
            imageYCbCr.at(y).at(x) = {
                Y_inverse.at(y).at(x),
                Cb_scaled_inverse.at(y/2).at(x/2),
                Cr_scaled_inverse.at(y/2).at(x/2)
            };
        }
    }
    input_stream.flush_to_byte();
    input_file.close();

    bitmap_image output_image {width,height};

    for (unsigned int y = 0; y < height; y++){
        for (unsigned int x = 0; x < width; x++){
            auto pixel_rgb = imageYCbCr.at(y).at(x).to_rgb();
            auto [r,g,b] = pixel_rgb;
            output_image.set_pixel(x,y,r,g,b);
        }
    }
    
    output_image.save_image(output_filename);
    
    return 0;
}