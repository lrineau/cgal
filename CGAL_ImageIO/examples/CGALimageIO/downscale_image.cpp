#include <CGAL/ImageIO.h>
#include <iostream>

#define SHOW(attribut) "\n  "#attribut": " << image->attribut
#define SHOWENUM(enumitem) #enumitem"=" << enumitem

int main(int argc, char** argv)
{
  if(argc == 0) return 1;

  _image* image = ::_readImage(argv[1]);
  if(image)
  {
    std::cerr
      << "Image infos:"
      << "\ndimensions"
      << SHOW(xdim)
      << SHOW(ydim)
      << SHOW(zdim)
      << SHOW(vdim)
      << "\nvoxel size"
      << SHOW(vx)
      << SHOW(vy)
      << SHOW(vz)
      << "\nword size (in bytes)"
      << SHOW(wdim)
      << "\n";
  }
  else
    std::cerr << "\"" << argv[1] << "\" is not a supported file.\n";

  if(argc != 3) return 1;
  const int xdim = image->xdim;
  const int ydim = image->ydim;
  const int zdim = image->zdim;
  const int wdim = image->wdim;

  _image image2 = *image;
  image2.xdim = xdim / 2;
  image2.ydim = ydim / 2;
  image2.zdim = zdim / 2;
  image2.vx = 2 * image->vx;
  image2.vy = 2 * image->vy;
  image2.vz = 2 * image->vz;
  image2.data =  ImageIO_alloc(std::size_t(wdim)*image2.xdim*image2.ydim*image2.zdim);
  const char* const ptr = (const char*)(image->data);
  char* out_ptr =               (char*)(image2.data);
  for(int k = 0; k < image2.zdim; ++k) {
    for(int j = 0; j < image2.ydim; ++j) {
      for(int i = 0; i < image2.xdim; ++i) {
        const char * const tmp_ptr = ptr +
          wdim * (2*i + (xdim * (2*j + ydim * std::size_t(2*k))));
        out_ptr = std::copy(tmp_ptr, tmp_ptr+wdim, out_ptr);
      }
    }
  }
  ::_writeImage(&image2, argv[2]);
  ImageIO_free(image2.data);
  ::_freeImage(image);
  return 0;
}
