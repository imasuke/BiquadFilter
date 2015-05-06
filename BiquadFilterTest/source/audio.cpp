#if _WIN32
#include "audio.hpp"
#else
#include "audio.hxx"
#endif
#include <fstream>

const int Audio::MAX_VALUE = 32768;

Audio::Audio(void)
: SignalT(0, 0, 0)
{
}
Audio::Audio(std::size_t spFreq, std::size_t qBit, std::size_t length)
: SignalT(spFreq, qBit, length)
{
}
Audio::Audio(const Audio &audio)
{
  spFreq_ = audio.spFreq_;
  qBit_ = audio.qBit_;
  length_ = audio.length_;
  samples_ = audio.samples_;
  filename = audio.filename;
}
Audio &Audio::operator=(const Audio &audio)
{
  spFreq_ = audio.spFreq_;
  qBit_ = audio.qBit_;
  length_ = audio.length_;
  samples_ = audio.samples_;
  filename = audio.filename;
  return (*this);
}
Audio::~Audio(void)
{
}

int Audio::open(const std::string &filename)
{
  short data;
  WaveFileHeader header;
  std::ifstream inputAudio;
  inputAudio.open(filename, std::ios::in | std::ios::binary);
  if (!inputAudio) {
    std::cerr << "audio file not exist" << std::endl;
    //throw std::exception("throw");
    return 0;
  }
  // set file name
  this->filename = filename;
  // read file header
  inputAudio.read(reinterpret_cast<char*>(&header), sizeof(header));

  // header information
  spFreq_ = header.fmtchunk.rate;
  qBit_ = header.fmtchunk.bits;
  length_ = header.datachunk.bytes / 2;
  samples_.resize(length_);

  // error process
  if(qBit_ == 8){
	std::cerr << "this function handles only 16 bits wav file" << std::endl;
	return 0;
  }

  // read data(this code is too slow)
  for (std::size_t n = 0; n < samples_.size(); ++n) {
    inputAudio.read(reinterpret_cast<char*>(&data), sizeof(short));
    samples_[n].time_ = n * dt();
    //samples_[n].sample_value_ = data / static_cast<double>(MAX_VALUE);
	samples_[n].sample_value_ = data;
  }

  return 0;
}
int Audio::save(const std::string &filename)
{
  std::ofstream outputAudio;
  outputAudio.open(filename, std::ios::out | std::ios::binary);
  if (!outputAudio) {
    std::cout << "Exception at: " << __FUNCTION__ << std::endl;
    //throw std::exception("Reason: Cannot Open File.");
    return 0;
  }

  char riff_chunk_ID[4];
  long riff_chunk_size;
  char riff_form_type[4];
  char fmt_chunk_ID[4];
  long fmt_chunk_size;
  short fmt_wave_format_type;
  short fmt_channel;
  long fmt_samples_per_sec;
  long fmt_bytes_per_sec;
  short fmt_block_size;
  short fmt_quantizationBitRate_per_sample;
  char data_chunk_ID[4];
  long data_chunk_size;
  short data;
//  double s;

  riff_chunk_ID[0] = 'R';
  riff_chunk_ID[1] = 'I';
  riff_chunk_ID[2] = 'F';
  riff_chunk_ID[3] = 'F';

  riff_form_type[0] = 'W';
  riff_form_type[1] = 'A';
  riff_form_type[2] = 'V';
  riff_form_type[3] = 'E';

  fmt_chunk_ID[0] = 'f';
  fmt_chunk_ID[1] = 'm';
  fmt_chunk_ID[2] = 't';
  fmt_chunk_ID[3] = ' ';

  fmt_chunk_size = 16;
  fmt_wave_format_type = 1;

  fmt_samples_per_sec = spFreq_;
  fmt_quantizationBitRate_per_sample = qBit_;

  data_chunk_ID[0] = 'd';
  data_chunk_ID[1] = 'a';
  data_chunk_ID[2] = 't';
  data_chunk_ID[3] = 'a';

  riff_chunk_size = 36 + length_ * 2;
  fmt_channel = 1;
  fmt_bytes_per_sec = spFreq_ * qBit_ / 8;
  fmt_block_size = qBit_ / 8;
  data_chunk_size = length_ * 2;

  // write header
  outputAudio.write(riff_chunk_ID, 4);
  outputAudio.write(reinterpret_cast<char*>(&riff_chunk_size), 4);
  outputAudio.write(riff_form_type, 4);
  outputAudio.write(fmt_chunk_ID, 4);
  outputAudio.write(reinterpret_cast<char*>(&fmt_chunk_size), 4);
  outputAudio.write(reinterpret_cast<char*>(&fmt_wave_format_type), 2);
  outputAudio.write(reinterpret_cast<char*>(&fmt_channel), 2);
  outputAudio.write(reinterpret_cast<char*>(&fmt_samples_per_sec), 4);
  outputAudio.write(reinterpret_cast<char*>(&fmt_bytes_per_sec), 4);
  outputAudio.write(reinterpret_cast<char*>(&fmt_block_size), 2);
  outputAudio.write(reinterpret_cast<char*>(&fmt_quantizationBitRate_per_sample), 2);
  outputAudio.write(data_chunk_ID, 4);
  outputAudio.write(reinterpret_cast<char*>(&data_chunk_size), 4);

  // write data
  for (std::size_t n = 0; n<length_; ++n) {
//    s = (samples_[n].sample_value_ + 1.0) / 2.0 * (2 * MAX_VALUE);
//    if (s > 2 * MAX_VALUE - 1) s = MAX_VALUE - 1;
//    if (s < 0) s = 0;
//    data = static_cast<short>((s + 0.5) - MAX_VALUE);
	  data = static_cast<short>(samples_[n].sample_value_);
	  outputAudio.write(reinterpret_cast<char*>(&data), sizeof(short));
  }

  return 0;
}

namespace AudioStream {
  /*------------------------------*/
  // AudioStream decoder class
  /*------------------------------*/
  decoder::decoder(void)
    : spFreq_(0), qBit_(0), BlockAlign_(0),
    nChannels_(0), AvgBytesPerSec_(0), DataLength_(0),
    write_cursole_(0), buffer_cursole_(0), frame_length_(0), frame_shift_(0)
  {
  }
  decoder::~decoder(void)
  {
  }

  void decoder::open(const std::string &filename)
  {
    Audio::WaveFileHeader header;

    // ファイルストリームのオープン
    input_.open(filename, std::ios::in | std::ios::binary);

    // ヘッダ読み取り
    input_.read(reinterpret_cast<char*>(&header), sizeof(header));

    // ヘッダ情報の提示
    std::cout << "----- Information of this wave file -----" << std::endl;
    std::cout << "Format ID     : " << header.fmtchunk.formatID << std::endl;
    std::cout << "Channel       : " << header.fmtchunk.channels << std::endl;
    std::cout << "Sampling Rate : " << header.fmtchunk.rate << std::endl;
    std::cout << "Bytes / Sec   : " << header.fmtchunk.velocity << std::endl;
    std::cout << "Block Size    : " << header.fmtchunk.blocksize << std::endl;
    std::cout << "Bit / Sample  : " << header.fmtchunk.bits << std::endl;
    std::cout << "Raw Data Size : " << header.datachunk.bytes << std::endl;
    std::cout << "-----------------------------------------" << std::endl;

    // WAVE情報
    spFreq_ = header.fmtchunk.rate;
    qBit_ = header.fmtchunk.bits;
    nChannels_ = header.fmtchunk.channels;
    BlockAlign_ = qBit_ / 8 * nChannels_;
    AvgBytesPerSec_ = spFreq_ * BlockAlign_;
    BufferLength_ = AvgBytesPerSec_ / 2;
    DataLength_ = header.datachunk.bytes / 2;

    // WAVE情報の提示
    std::cout << "SamplePerSec   = " << spFreq_ << std::endl;
    std::cout << "BitsPerSample  = " << qBit_ << std::endl;
    std::cout << "BlockAlign     = " << BlockAlign_ << std::endl;
    std::cout << "AvgBytesPerSec = " << AvgBytesPerSec_ << std::endl;
    std::cout << "BufferLength   = " << BufferLength_ << std::endl;
    std::cout << "DataLength     = " << DataLength_ << std::endl;

    // バッファリサイズ（1秒）
    bufferPerSec_.first.resize(BufferLength_);
    bufferPerSec_.second.resize(BufferLength_);

    // 初期位置(ヘッダ)
    write_cursole_ += sizeof(Audio::WaveFileHeader);

    return;
  }

  void decoder::stockBuffer(void)
  {
    beforeBuffer_ = bufferPerSec_.first;
  }
  void decoder::addNextBuffer(void)
  {
    std::vector<short> &data = bufferPerSec_.second;
    int write_cursole = (static_cast<std::size_t>(input_.tellg())
      - sizeof(Audio::WaveFileHeader)) / sizeof(short);
    if (write_cursole + BufferLength_ < DataLength_) {
      for (std::size_t n = 0; n < BufferLength_; ++n) {
        input_.read(reinterpret_cast<char*>(&data[n]), sizeof(short));
        write_cursole_ += sizeof(short);
      }
    }
    else {
      std::size_t len = DataLength_ - write_cursole;
      for (std::size_t n = 0; n < len; ++n) {
        input_.read(reinterpret_cast<char*>(&data[n]), sizeof(short));
        write_cursole += sizeof(short);
      }
    }
    buffer_cursole_++;
    flipBuffer();
  }
  void decoder::flipBuffer(void)
  {
    std::swap(bufferPerSec_.first, bufferPerSec_.second);
  }

  void decoder::seek(std::size_t n)
  {
    int write_cursole = (static_cast<std::size_t>(input_.tellg()) - sizeof(Audio::WaveFileHeader)) / sizeof(short);
    std::cout << write_cursole << "\t" << write_cursole_ / 2 << std::endl;
    input_.seekg(sizeof(Audio::WaveFileHeader) + sizeof(short)* n);
  }

  void decoder::setFrameLength(std::size_t len)
  {
    frame_length_ = len;
  }
  void decoder::setFrameShift(std::size_t shift)
  {
    frame_shift_ = shift;
  }

  double decoder::at(std::size_t n)
  {
    std::size_t block = n / BufferLength_;
    std::size_t index = n % BufferLength_;

    // read next buffer
    if (block >= buffer_cursole_) {
      do {
        stockBuffer();
        addNextBuffer();
      } while (block == buffer_cursole_);
    }
    double res;
    res = (bufferPerSec_.first[index]);
    if (buffer_cursole_ - 1 != block) {
      res = beforeBuffer_[index];
    }
    return res / static_cast<double>(Audio::MAX_VALUE);
  }
  double decoder::operator()(std::size_t m, std::size_t t)
  {
    return at(m*frame_shift_+t);
  }

  void decoder::extract(Audio *out, std::size_t sample1, std::size_t sample2)
  {
    if (sample1 > sample2) std::swap(sample1, sample2);
    std::size_t len = sample2 - sample1;
    for (std::size_t n = 0; n < len; ++n) {
      (*out)[n] = this->at(n + sample1);
    }
  }
  void decoder::frame(Audio *out, std::size_t nframe)
  {
    for (std::size_t n = 0; n < frame_length_; ++n) {
      (*out)[n] = this->at(n + nframe*frame_shift_);
    }
  }

  /*------------------------------*/
  // AudioStream encoder class
  /*------------------------------*/
  encoder::encoder(void)
    : spFreq_(0), qBit_(0), BlockAlign_(0),
    nChannels_(0), AvgBytesPerSec_(0), DataLength_(0)
  {
  }
  encoder::~encoder(void)
  {
  }
  int encoder::open(const std::string &filename,
    std::size_t DataLength,
    std::size_t SamplePerSec,
    std::size_t BitsPerSample)
  {
    // audio information
    spFreq_ = SamplePerSec;
    qBit_ = BitsPerSample;
    DataLength_ = DataLength;
    nChannels_ = 1;
    BlockAlign_ = qBit_ / 8 * nChannels_;
    AvgBytesPerSec_ = spFreq_ * BlockAlign_;
    BufferLength_ = AvgBytesPerSec_ / 2;

    // open file stream
    outputAudio_.open(filename, std::ios::out | std::ios::binary);
    if (!outputAudio_) {
      std::cout << "Exception at: " << __FUNCTION__ << std::endl;
      //throw std::exception("Reason: cannot open file");
      return -1;
    }

    // header information
    char riff_chunk_ID[4];
    long riff_chunk_size;
    char riff_form_type[4];
    char fmt_chunk_ID[4];
    long fmt_chunk_size;
    short fmt_wave_format_type;
    short fmt_channel;
    long fmt_samples_per_sec;
    long fmt_bytes_per_sec;
    short fmt_block_size;
    short fmt_quantizationBitRate_per_sample;
    char data_chunk_ID[4];
    long data_chunk_size;

    riff_chunk_ID[0] = 'R';
    riff_chunk_ID[1] = 'I';
    riff_chunk_ID[2] = 'F';
    riff_chunk_ID[3] = 'F';

    riff_form_type[0] = 'W';
    riff_form_type[1] = 'A';
    riff_form_type[2] = 'V';
    riff_form_type[3] = 'E';

    fmt_chunk_ID[0] = 'f';
    fmt_chunk_ID[1] = 'm';
    fmt_chunk_ID[2] = 't';
    fmt_chunk_ID[3] = ' ';

    fmt_chunk_size = 16;
    fmt_wave_format_type = 1;

    fmt_samples_per_sec = spFreq_;  // sampling frequency
    fmt_quantizationBitRate_per_sample = qBit_; // quantinzation bit

    data_chunk_ID[0] = 'd';
    data_chunk_ID[1] = 'a';
    data_chunk_ID[2] = 't';
    data_chunk_ID[3] = 'a';

    riff_chunk_size = 36 + DataLength_ * 2;
    fmt_channel = 1;
    fmt_bytes_per_sec = spFreq_ * qBit_ / 8;
    fmt_block_size = qBit_ / 8;
    data_chunk_size = DataLength_ * 2;

    outputAudio_.write(riff_chunk_ID, 4);
    outputAudio_.write(reinterpret_cast<char*>(&riff_chunk_size), 4);
    outputAudio_.write(riff_form_type, 4);
    outputAudio_.write(fmt_chunk_ID, 4);
    outputAudio_.write(reinterpret_cast<char*>(&fmt_chunk_size), 4);
    outputAudio_.write(reinterpret_cast<char*>(&fmt_wave_format_type), 2);
    outputAudio_.write(reinterpret_cast<char*>(&fmt_channel), 2);	//
    outputAudio_.write(reinterpret_cast<char*>(&fmt_samples_per_sec), 4);
    outputAudio_.write(reinterpret_cast<char*>(&fmt_bytes_per_sec), 4);
    outputAudio_.write(reinterpret_cast<char*>(&fmt_block_size), 2);
    outputAudio_.write(reinterpret_cast<char*>(&fmt_quantizationBitRate_per_sample), 2);
    outputAudio_.write(data_chunk_ID, 4);
    outputAudio_.write(reinterpret_cast<char*>(&data_chunk_size), 4);

    return 0;
  }
  void encoder::at(std::size_t n, double value)
  {
    short data;
    double s;
    s = (value + 1.0) / 2.0 * (2*Audio::MAX_VALUE);
    if (s > (2 * Audio::MAX_VALUE)) s = (2 * Audio::MAX_VALUE);
    if (s < 0.0) s = 0.0;
    data = static_cast<short>((s + 0.5) - Audio::MAX_VALUE);
    outputAudio_.write(reinterpret_cast<char*>(&data), sizeof(short));
  }
}
