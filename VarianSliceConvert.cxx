#include "itkImage.h"
#include "itkImageIORegion.h"
#include <iostream>

#include "itkByteSwapper.h"
#include "itkRGBPixel.h"
#include "itkRGBAPixel.h"
#include "vnl/vnl_matrix.h"
#include <cstdio>
#include <fstream>


//#include "itkFDF2CommonImageIO.h"
//#include "itkFDF2ImageIOFactory.h"
//#include "itkFDF2ImageIO.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

// Remove a particular type of character from a string
std::string
RemoveCharacters(std::string line, char character)
{
  line.erase(std::remove(line.begin(), line.end(), character), line.end());
  return line;
}

void
Tokenize(const std::string & str, std::vector<std::string> & tokens, const std::string & delimiters)
{
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  std::string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos)
  {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}

std::string
ParseLine(std::string line)
{
  // strip *
  line = RemoveCharacters(line, '*');
  line = RemoveCharacters(line, '\"');
  line = RemoveCharacters(line, '[');
  line = RemoveCharacters(line, ']');
  line = RemoveCharacters(line, '\r');

  // Need to deal with space between {}
  std::string::size_type startBracketPosition = line.find_first_of("{", 0);
  std::string::size_type endBracketPosition = line.find_first_of("}", startBracketPosition);

  if (startBracketPosition != std::string::npos && endBracketPosition != std::string::npos)
  {
    std::string element = line.substr(startBracketPosition, endBracketPosition - startBracketPosition);

    // Find whitespace within {} and erase
    std::string::size_type whiteSpacePosition = line.find_first_of(" ", startBracketPosition);

    while (whiteSpacePosition != std::string::npos)
    {
      line.erase(whiteSpacePosition, 1);
      whiteSpacePosition = line.find_first_of(" ", whiteSpacePosition);
    }
  }

  return line;
}

template <typename T>
void
ConvertFromString(std::string s, T & value)
{
  std::stringstream str;
  str << s;
  str >> value;
}

template <typename T>
void
StringToVector(std::string value, std::vector<T> & values)
{
  std::vector<std::string> tokens;

  // value consists of something like {256,256}
  std::string::size_type startBracketPosition = value.find_first_of("{", 0);
  std::string::size_type endBracketPosition = value.find_first_of("}", startBracketPosition);

  if (startBracketPosition != std::string::npos && endBracketPosition != std::string::npos)
  {
    std::string elements = value.substr(startBracketPosition + 1, endBracketPosition - startBracketPosition - 1);


    Tokenize(elements, tokens, ",");
  }

  T element;

  for (auto & token : tokens)
  {
    ConvertFromString(token, element);
    values.push_back(element);
  }
}

template <typename T>
void
PrintVector(std::ostream & os, std::string name, const std::vector<T> & vect)
{
  int size = vect.size();

  os << name << " {";

  for (int i = 0; i < size; i++)
  {
    os << vect[i];

    if (i < size - 1)
      os << ", ";
  }

  os << "}" << std::endl;
}


bool
CanReadFile(const std::string filename)
{

  if (filename.empty())
  {
    return false;
  }

  bool                   extensionFound = false;
  std::string::size_type FDF2Pos = filename.rfind(".FDF");
  if ((FDF2Pos != std::string::npos) && (FDF2Pos == filename.length() - 4))
  {
    extensionFound = true;
  }

  FDF2Pos = filename.rfind(".fdf");
  if ((FDF2Pos != std::string::npos) && (FDF2Pos == filename.length() - 4))
  {
    extensionFound = true;
  }

  if (!extensionFound)
  {
    //itkDebugMacro(<< "The filename extension is not recognized");
    return false;
  }

  std::ifstream inFile;
  inFile.open(filename.c_str(), std::ios::in | std::ios::binary);
  if (!inFile)
  {
    return false;
  }

  // Check for a neccessary header variable

  inFile.close();
  return true;
}

std::string
GetStorageType( std::string filename )
{
  std::string storage;

  if (!CanReadFile(filename)) {
    return storage;
  }

  std::string              line;
  std::vector<std::string> tokens;
  std::string              type, name, value;

  std::ifstream inFile(filename.c_str(), std::ios::in | std::ios::binary);

  // Check if there was an error opening the file
  if (!inFile)
  {
    std::cout << "Unable to open the file\n";
    return storage;
  }

  while (getline(inFile, line, '\n'))
  {
    if (line == "\0")
    {
      break;
    }

    if (line.find('\f') != std::string::npos)
    {
      //std::cout << "Line has form feed" << std::endl;
      char b;
      inFile.read(&b,1);
      while (b != 0) {
        inFile.read(&b,1);
      }
      break;
    }

    // Formats the lines in the FDF header such as removing whitespace between {}
    line = ParseLine(line);
    Tokenize(line, tokens, " ;");

    if (tokens.size() >= 4)
    {
      //std::cout << "tokens[0]: " << tokens[0] << std::endl;
      type = tokens[0];
      name = tokens[1];
      value = tokens[3];

      // Get the binary data type
      if (name == "storage")
      {
        return value;
      }
    }

    tokens.clear();

  }

  inFile.close();
  return storage;
}
unsigned int
GetComponentSize(std::string storage)
{
  if (storage=="unsigned char")
      return sizeof(unsigned char);
  if (storage=="char")
      return sizeof(char);
  if (storage=="unsigned short")
      return sizeof(unsigned short);
  if (storage=="short")
      return sizeof(short);
  if (storage=="unsigned int")
      return sizeof(unsigned int);
  if (storage=="int")
      return sizeof(int);
  if (storage=="unsigned long")
      return sizeof(unsigned long);
  if (storage=="long")
      return sizeof(long);
  if (storage=="unsigned long long")
      return sizeof(unsigned long long);
  if (storage=="long long")
      return sizeof(long long);
  if (storage=="float")
      return sizeof(float);
  if (storage=="double")
      return sizeof(double);

  std::cerr << "Unknown component type: " << storage << std::endl;
  return 0;
}

void
SwapBytesIfNecessary(void * buffer, unsigned long numberOfPixels, std::string byteOrder, std::string storage)
{

  if (storage=="char")
  {
    if (byteOrder == "littleendian")
    {
      itk::ByteSwapper<char>::SwapRangeFromSystemToLittleEndian((char *)buffer, numberOfPixels);
    }
    else if (byteOrder == "bigendian")
    {
      itk::ByteSwapper<char>::SwapRangeFromSystemToBigEndian((char *)buffer, numberOfPixels);
    }
    return;
  }
  if (storage=="float")
  {
    if (byteOrder == "littleendian")
    {
      itk::ByteSwapper<float>::SwapRangeFromSystemToLittleEndian((float *)buffer, numberOfPixels);
    }
    else if (byteOrder == "bigendian")
    {
      itk::ByteSwapper<float>::SwapRangeFromSystemToBigEndian((float *)buffer, numberOfPixels);
    }
    return;
  }
  if (storage=="unsigned char")
  {
    if (byteOrder == "littleendian")
    {
      itk::ByteSwapper<unsigned char>::SwapRangeFromSystemToLittleEndian((unsigned char *)buffer, numberOfPixels);
    }
    else if (byteOrder == "bigendian")
    {
      itk::ByteSwapper<unsigned char>::SwapRangeFromSystemToBigEndian((unsigned char *)buffer, numberOfPixels);
    }
    return;
  }
  if ( storage=="short")
  {
    if (byteOrder == "littleendian")
    {
      itk::ByteSwapper<short>::SwapRangeFromSystemToLittleEndian((short *)buffer, numberOfPixels);
    }
    else if (byteOrder == "bigendian")
    {
      itk::ByteSwapper<short>::SwapRangeFromSystemToBigEndian((short *)buffer, numberOfPixels);
    }
    return;
  }
  if (storage=="unsigned short")
  {
    if (byteOrder == "littleendian")
    {
      itk::ByteSwapper<unsigned short>::SwapRangeFromSystemToLittleEndian((unsigned short *)buffer, numberOfPixels);
    }
    else if (byteOrder == "bigendian")
    {
      itk::ByteSwapper<unsigned short>::SwapRangeFromSystemToBigEndian((unsigned short *)buffer, numberOfPixels);
    }
    return;
  }
  if (storage=="int")
  {
    if (byteOrder == "littleendian")
    {
      itk::ByteSwapper<int>::SwapRangeFromSystemToLittleEndian((int *)buffer, numberOfPixels);
    }
    else if (byteOrder == "bigendian")
    {
      itk::ByteSwapper<int>::SwapRangeFromSystemToBigEndian((int *)buffer, numberOfPixels);
    }
    return;
  }
  if (storage=="unsigned int")
  {
    if (byteOrder == "littleendian")
    {
      itk::ByteSwapper<unsigned int>::SwapRangeFromSystemToLittleEndian((unsigned int *)buffer, numberOfPixels);
    }
    else if (byteOrder == "bigendian")
    {
      itk::ByteSwapper<unsigned int>::SwapRangeFromSystemToBigEndian((unsigned int *)buffer, numberOfPixels);
    }
    return;
  }
  if (storage=="long")
  {
    if (byteOrder == "littleendian")
    {
      itk::ByteSwapper<long>::SwapRangeFromSystemToLittleEndian((long *)buffer, numberOfPixels);
    }
    else if (byteOrder == "bigendian")
    {
      itk::ByteSwapper<long>::SwapRangeFromSystemToBigEndian((long *)buffer, numberOfPixels);
    }
    return;
  }
  if (storage=="unsigned long")
  {
    if (byteOrder == "littleendian")
    {
      itk::ByteSwapper<unsigned long>::SwapRangeFromSystemToLittleEndian((unsigned long *)buffer, numberOfPixels);
    }
    else if (byteOrder == "bigendian")
    {
      itk::ByteSwapper<unsigned long>::SwapRangeFromSystemToBigEndian((unsigned long *)buffer, numberOfPixels);
    }
    return;
  }

}


template< typename ImageType >
bool
ReadImage( std::string filename, typename ImageType::Pointer img )
{

  if (!CanReadFile(filename)) {
    return false;
  }

  const unsigned numDim = ImageType::ImageDimension;

  std::string              line;
  std::vector<std::string> tokens;
  std::string              type, name, value;

  typename ImageType::RegionType region;
  std::vector<double> roi;
  std::string storage;
  std::string byteOrder;

  std::ifstream inFile(filename.c_str(), std::ios::in | std::ios::binary);

  // Check if there was an error opening the file
  if (!inFile)
  {
    std::cout << "Unable to open the file\n";
    return false;
  }

  while (getline(inFile, line, '\n'))
  {
    if (line == "\0")
    {
      break;
    }

    if (line.find('\f') != std::string::npos)
    {
      //std::cout << "Line has form feed" << std::endl;
      char b;
      inFile.read(&b,1);
      while (b != 0) {
        inFile.read(&b,1);
      }
      break;
    }

    // Formats the lines in the FDF header such as removing whitespace between {}
    line = ParseLine(line);
    Tokenize(line, tokens, " ;");

    if (tokens.size() >= 4)
    {
      type = tokens[0];
      name = tokens[1];
      value = tokens[3];

      if (name == "spatial_rank")
      {
        //this->m_SpatialRank = value;
      }

      if (name == "matrix")
      {
        std::vector<float> dimensions;
        StringToVector(value, dimensions);

        typename ImageType::RegionType::SizeType  size;
        typename ImageType::RegionType::IndexType index;

        for (unsigned int i = 0; i < numDim; i++)
        {
          size[i] = dimensions[i];
          index[i] = 0;
        }

        region.SetSize(size);
        region.SetIndex(index);
        img->SetRegions(region);
      }

      if (name == "orientation")
      {
        std::vector<double> orientation;
        StringToVector(value, orientation);

        vnl_matrix<double> testDirections(numDim, numDim);
        typename ImageType::DirectionType dirMat;

        for (unsigned int i = 0; i < numDim; i++)
        {
          std::vector<double> componentVector;
          for (unsigned int j = 0; j < numDim; j++)
          {
            double val = orientation[i * numDim + j];
            testDirections(j, i) = val;
            dirMat(j,i) = val;
            componentVector.push_back(val);
          }

        }
        // check for degenerate dimensions. this will happen
        // if the dimension of the image is 2 but the
        // direction matrix in the file is 3x3.
        // if direction matrix is degenerate, punt and set
        // directions to identity
        if (vnl_determinant(testDirections) == 0)
        {
          for (unsigned int i = 0; i < numDim; i++)
          {
            std::vector<double> componentVector;
            for (unsigned int j = 0; j < numDim; j++)
            {
              double val = i == j ? 1.0 : 0.0;
              componentVector.push_back(val);
              dirMat(i,j) = val;
            }
          }
        }


        img->SetDirection(dirMat);
      }
      if (name == "span")
      {
        std::vector<double> span;
        StringToVector(value, span);
      }

      if (name == "origin")
      {
        std::vector<float> origin;
        StringToVector(value, origin);

        //if (this->GetNumberOfDimensions() < origin.size())
        //{
        //img->SetNumberOfDimensions(origin.size());
        //}

        typename ImageType::PointType itkorigin;
        for (unsigned int i = 0; i < origin.size(); i++)
        {
          itkorigin[i] =  origin[i] / 10.0;
        }
        img->SetOrigin(itkorigin);
      }

      if (name == "roi")
      {
        StringToVector(value, roi);
      }

      if (name == "location")
      {
        std::vector<double> location;
        StringToVector(value, location);
      }

      if (name == "bigendian")
      {
        if (value == "0")
        {
          byteOrder = "littleendian";
        }
        else
        {
          byteOrder = "bigendian";
        }
      }

      // Get the binary data type
      if (name == "storage")
      {
        storage = value;
      }

      // Get the bits
      if (name == "bits")
      {
        //std::vector<double> bits;
        //ConvertFromString(value, bits);
      }

      // Get the checksum
      if (name == "checksum")
      {
        //std::vector<double> checksum;
        //ConvertFromString(value, checksum);
      }

    }

    tokens.clear();
  }

  typename ImageType::SpacingType spacing;
  for (unsigned int i = 0; i < numDim; i++)
  {
    spacing[i] = (roi[i] * 10) / region.GetSize()[i];
  }

  inFile.seekg(0, std::ios::end);
  long int fileSize = inFile.tellg();
  unsigned long sizeInPixels = 1;
  for ( unsigned i=0; i<numDim; i++) {
    sizeInPixels *= region.GetSize()[i];
  }
  unsigned long sizeInBytes = sizeInPixels *= GetComponentSize( storage );
  inFile.seekg(fileSize - sizeInBytes);
  std::cout << "Size in bytes: " << sizeInBytes << std::endl;

  img->Allocate();
  //unsigned char * buffer = img->GetBufferPointer();
  auto * p = reinterpret_cast<char *>(img->GetBufferPointer());
  inFile.read(p, sizeInBytes);

  bool success = !inFile.bad();
  inFile.close();
  if (!success)
  {
    std::cerr << "Error reading image data." << std::endl;
  }



  void * buffer = reinterpret_cast<void *>(img->GetBufferPointer());
  SwapBytesIfNecessary(buffer, sizeInPixels, byteOrder, storage);


  return true;
}


int main(int argc, char * argv[])
{

  if ( argc != 3 ) {
    std::cout << "Usage: VarianConvert input.fdf output.nii.gz" << std::endl;
  }


  //using ReaderType = itk::ImageFileReader<ImageType>;

  //itk::FDF2ImageIOFactory::RegisterOneFactory();


  std::string filename(argv[1]);
  if ( CanReadFile( filename ) ) {
    std::cout << "Can read the input file: " << filename << std::endl;
  }
  else {
    std::cout << "Unable to read the input file: " << filename << std::endl;
  }

  std::string storage = GetStorageType( filename );
  std::cout << "ValueType: " << storage << std::endl;

  if ( storage == "float" )
  {
    using ImageType = itk::Image< float, 2 >;
    using WriterType = itk::ImageFileWriter<ImageType>;
    using ImagePointer = ImageType::Pointer;
    ImagePointer img = ImageType::New();

    if ( ReadImage<ImageType>( filename, img ) ) {
      std::cout << "Can read header: " << filename << std::endl;
      std::cout << img << std::endl;

      typename WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( argv[2] );
      writer->SetInput(img);
      writer->Update();

    }
    else {
      std::cout << "Unable to read header: " << filename << std::endl;
    }

  }





  return 0;


  //std::string filename(argv[1]);

  /*
  std::ifstream inFile1;
  inFile1.open(filename, std::ios::in | std::ios::binary);
  if (!inFile1)
  {
    std::string  msg = "File \"" + filename + "\" cannot be opened.";
    std::cout << msg << std::endl;
  }


  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  bool success=false;
  try {
    reader->Update();
    success=true;
  }
  catch( itk::ExceptionObject & err ) {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
  }

  if (success)
  {
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( argv[2] );
    writer->SetInput( reader->GetOutput() );
    writer->Update();
    return 0;
  }

  return 0;
  */
}
