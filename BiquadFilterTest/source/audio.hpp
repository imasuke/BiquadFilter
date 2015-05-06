#ifndef _AUDIO_HEADER_
#define _AUDIO_HEADER_

#if _WIN32
#include "signalT.hpp"
#else
#include "signalT.hxx"
#endif
#include <fstream>

class Audio : public SignalT < double >
{
public:
	static const int MAX_VALUE; // max = 32768

	// WAVE format chunk
	typedef struct {
		unsigned int   id;       // "fmt "
		unsigned int   bytes;    // fmt chunk bytes
		unsigned short formatID; // format ID
		unsigned short channels; // channel number
		unsigned int   rate;     // sampling rate
		unsigned int   velocity; // data velocity
		unsigned short blocksize;// block size
		unsigned short bits;     // quantization bit number
	} fmtChunk;
	// WAVE data chunk
	typedef struct {
		unsigned int id;       // "data"
		unsigned int bytes;    // signal data bytes
	} dataChunk;
	// WAVE file(RIFF) header
	typedef struct {
		unsigned int riffID; // "riff"
		unsigned int size_8; // file size-8
		unsigned int wavID;  // "WAVE"
		fmtChunk fmtchunk;
		dataChunk datachunk;
	} WaveFileHeader;

	std::string filename;

public:

	// constructor
	Audio(void);
	Audio(std::size_t spFreq, std::size_t qBit, std::size_t length);
	Audio(const std::string &filename);
	Audio(const Audio &audio);
	Audio &operator=(const Audio &audio);
	~Audio(void);

	// i/o method
	int open(const std::string &filename);
	int save(const std::string &filename);

};


namespace AudioStream {

	typedef std::pair< std::vector<short>, std::vector<short> > Buffer;

	class decoder {
	private:
		std::size_t spFreq_;    // sampling frequency
		std::size_t qBit_;      // quantized bits

		std::size_t BlockAlign_;
		std::size_t nChannels_;
		std::size_t AvgBytesPerSec_;  // don't need?
		std::size_t BufferLength_;    // don't need?
		std::size_t DataLength_;

		// data seek cursole
		std::size_t write_cursole_;
		std::size_t buffer_cursole_;

		// buffer
		std::vector<short> beforeBuffer_;
		Buffer bufferPerSec_;
		std::ifstream input_;

		std::size_t frame_length_;
		std::size_t frame_shift_;

	public:
		decoder(void);
		~decoder(void);

		// i/o(decoder class only open audio)
		void open(const std::string &filename);

		// buffer
		void stockBuffer(void);   // stock first buffer
		void flipBuffer(void);    // flip primary and secondary buffer
		void addNextBuffer(void); // add next buffer

		void seek(std::size_t n);

		// set parameter
		void setFrameLength(std::size_t len);
		void setFrameShift(std::size_t shift);

		// get wave parameter
		std::size_t spFreq(void)
		{
			return spFreq_;
		}
		std::size_t qBit(void)
		{
			return qBit_;
		}
		std::size_t length(void)
		{
			return DataLength_;
		}

		// accesser
		double at(std::size_t n);
		double operator()(std::size_t m, std::size_t t);

		// extract
		void extract(Audio *out, std::size_t sample1, std::size_t sample2);
		void frame(Audio *out, std::size_t nframe);

		// utility method(to transfer audio lib?)
		double dt(void)
		{
			return 1.0 / static_cast<double>(spFreq_);
		}
		std::size_t number_of_frames(void)
		{
			return DataLength_ / frame_shift_;
		}
		std::size_t last_frame_size(void)
		{
			return DataLength_ - (number_of_frames() * frame_shift_);
		}
		std::size_t frame_length(void)
		{
			return frame_length_;
		}
		std::size_t frame_shift(void)
		{
			return frame_shift_;
		}
	};

	class encoder {
	private:
		std::size_t spFreq_;
		std::size_t qBit_;
		std::size_t BlockAlign_;
		std::size_t nChannels_;
		std::size_t AvgBytesPerSec_;
		std::size_t BufferLength_;
		std::size_t DataLength_;

		std::ofstream outputAudio_;

	public:
		encoder(void);
		~encoder(void);
		int open(const std::string &filename,
			std::size_t DataLength,
			std::size_t SamplePerSec,
			std::size_t BitsPerSample);
		void at(std::size_t n, double value);
	};
}

// pre-version AudioStream
#if 0
namespace AudioStream_pre {

	static const int FRAME_LENGTH = 512;	// frame length
	static const int FRAME_CYCLE = 128;	  // frame cycle

	typedef std::pair< std::vector<short>, std::vector<short> > Buffer;

	class decoder {
	private:

		std::size_t SamplePerSec_;
		std::size_t BitsPerSample_;
		std::size_t BlockAlign_;
		std::size_t nChannels_;
		std::size_t AvgBytesPerSec_;
		std::size_t BufferLength_;
		std::size_t DataLength_;

		std::size_t write_cursole_;
		std::size_t buffer_cursole_;


		std::vector<short> beforeBuffer_;
		Buffer bufferPerSec_;
		std::ifstream inputStream_;

		std::size_t frame_length_;
		std::size_t frame_cycle_;

	public:
		decoder(void)
			: SamplePerSec_(0), BitsPerSample_(0), BlockAlign_(0),
			nChannels_(0), AvgBytesPerSec_(0), DataLength_(0),
			write_cursole_(0), buffer_cursole_(0), frame_length_(0), frame_cycle_(0)
		{
		}
		~decoder(void)
		{
		}

		//// temp ----------------------------
		double SampleToTime(int sample)   { return static_cast<double>(sample) / static_cast<double>(SamplePerSec_); }			// sample point to time[sec]
		int    TimeToSample(double _time) { return static_cast<int>(_time * static_cast<double>(SamplePerSec_)); }			// time[sec] to sample point
		Audio extractTime(double time1, double time2)
		{
			int n1 = TimeToSample(time1), n2 = TimeToSample(time2);
			std::cout << n2 << ", " << n1 << std::endl;
			Audio res(SamplePerSec_, BitsPerSample_, n2 - n1);
			for (int n = n1; n < n2; n++) {
				res[n - n1] = (*this)[n];
			}
			return res;
		}
		////////////////////////////////

		void open(const std::string &filename)
		{
			Audio::WaveFileHeader header;

			// ファイルストリームのオープン
			inputStream_.open(filename, std::ios::in | std::ios::binary);

			// ヘッダ読み取り
			inputStream_.read(reinterpret_cast<char*>(&header), sizeof(header));

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
			SamplePerSec_ = header.fmtchunk.rate;
			BitsPerSample_ = header.fmtchunk.bits;
			nChannels_ = header.fmtchunk.channels;
			BlockAlign_ = BitsPerSample_ / 8 * nChannels_;
			AvgBytesPerSec_ = SamplePerSec_ * BlockAlign_;
			BufferLength_ = AvgBytesPerSec_ / 2;
			DataLength_ = header.datachunk.bytes / 2;

			// WAVE情報の提示
			std::cout << "SamplePerSec   = " << SamplePerSec_ << std::endl;
			std::cout << "BitsPerSample  = " << BitsPerSample_ << std::endl;
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

		std::size_t getSamplePerSec(void)  const { return SamplePerSec_; }
		std::size_t getBitsPerSample(void) const { return BitsPerSample_; }
		std::size_t getDataLength(void)    const { return DataLength_; }
		std::size_t getBufferLength(void)  const { return BufferLength_; }

		Audio ExtractAudio(int sample1, int sample2)
		{
			int length = sample2 - sample1;
			Audio res(SamplePerSec_, BitsPerSample_, sample2 - sample1);
			for (int n = 0; n < length; n++) {
				//std::cout << (*this)[n+sample1] << std::endl;
				res[n] = (*this)[n + sample1];
				//std::cin.get();	
			}
			return res;
		}

		// フレーム長の設定
		void setFrameLength(std::size_t size = FRAME_LENGTH)
		{
			frame_length_ = size;
		}

		// フレーム周期の設定
		void setFrameCycle(std::size_t cycle = FRAME_CYCLE)
		{
			frame_cycle_ = cycle;
		}

		// フレーム数の取得
		std::size_t NumberOfFrame(void) const
		{
			if (frame_cycle_ == 0) {
				std::cout << "Exception at: " << __FUNCTION__ << std::endl;
				//throw std::exception("Reason: Frame Cycle is not set.");
				return 0;
			}
			return this->DataLength_ / this->frame_cycle_;
		}

		// 最終フレームのサイズを取得
		std::size_t LastFrameSize(void) const
		{
			return DataLength_ - (NumberOfFrame() * frame_cycle_);
		}

		// フレーム長の取得
		std::size_t FrameLength(void)
		{
			return frame_length_;
		}

		// フレーム周期の取得
		std::size_t FrameCycle(void)
		{
			return frame_cycle_;
		}

		// 
		void stockBuffer(void)
		{
			beforeBuffer_ = bufferPerSec_.first;
		}

		// 次のバッファを読み込む
		void addNextBuffer(void)
		{
			// セカンダリバッファのエイリアス
			std::vector<short> &data = bufferPerSec_.second;

			// セカンダリバッファにデータを書き込み，カーソルを進める
			int write_cursole = (static_cast<std::size_t>(inputStream_.tellg()) - sizeof(Audio::WaveFileHeader)) / sizeof(short);

			if (write_cursole + BufferLength_ < DataLength_) {
				for (std::size_t n = 0; n < BufferLength_; ++n) {
					inputStream_.read((char*)&data[n], sizeof(short));
					write_cursole_ += sizeof(short);
				}
			}
			else {
				std::size_t len = DataLength_ - write_cursole;
				for (std::size_t n = 0; n < len; ++n) {
					inputStream_.read((char*)&data[n], sizeof(short));
					write_cursole_ += sizeof(short);
				}
			}

			// バッファ位置子の移動
			buffer_cursole_++;

			// プライマリバッファとセカンダリバッファの交換
			flipBuffer();

			return;
		}

		// プライマリとセカンダリを交換
		void flipBuffer(void)
		{
			std::swap(bufferPerSec_.first, bufferPerSec_.second);
		}

		// 指定サンプル値へシーク
		void seek(std::size_t n)
		{
			int write_cursole = (static_cast<std::size_t>(inputStream_.tellg()) - sizeof(Audio::WaveFileHeader)) / sizeof(short);
			std::cout << write_cursole << "\t" << write_cursole_ / 2 << std::endl;
			inputStream_.seekg(sizeof(Audio::WaveFileHeader) + sizeof(short)* n);
			std::cout << "tell: " << inputStream_.tellg() << std::endl;
			//std::cin.get();
		}

		// ランダムアクセス可能
		// !! 毎回シークするので読み込み速度は遅い
		double operator[](std::size_t n)
		{
			// 指定位置までシーク
			//int write_cursole = (static_cast<std::size_t>(inputStream_.tellg()) - sizeof(Audio::WaveFileHeader)) / sizeof(short);
			inputStream_.seekg(sizeof(Audio::WaveFileHeader) + sizeof(short)* n);

			short data;
			inputStream_.read(reinterpret_cast<char*>(&data), sizeof(short));

			return data / static_cast<double>(Audio::MAX_VALUE);
		}

		// 順次アクセスのみ対応
		// !! nを戻さないこと（windowingなどには不向き）
		double at(int n)
		{
			//	std::size_t index;
			//	index = n % BufferLength_;
			std::size_t block = n / BufferLength_;
			std::size_t index = n % BufferLength_;

			// 次のバッファを読み込む
			if (block >= buffer_cursole_) {
				do {
					stockBuffer();
					addNextBuffer();
				} while (block == buffer_cursole_);
			}

			double res;
			res = (bufferPerSec_.first[index]);
			if (buffer_cursole_ - 1 != block) res = beforeBuffer_[index];

			return res / static_cast<double>(Audio::MAX_VALUE);
		}

		// frame processing
		// s(m,t) : sample 't' in frame cycle 'm'
		double operator()(std::size_t m, std::size_t t)
		{
			//return (*this)[m*frame_cycle_+t];
			return at(m*frame_cycle_ + t);
		}
		void frame(std::size_t m, Audio *out)
		{
			if (out == NULL) {
				std::cout << "Exception at: " << __FUNCTION__ << std::endl;
				//throw std::exception("Reason: argument pointer is NULL");
				return;
			}

			out->allocate(SamplePerSec_, BitsPerSample_, frame_length_);

			for (std::size_t i = 0; i < frame_length_; ++i) {
				(*out)[i] = (*this)[i + m*frame_cycle_];
			}

		}

		double dt(void) const
		{
			return 1.0 / static_cast<double>(SamplePerSec_);
		}
	};

	class encoder {
	private:
		std::size_t SamplePerSec_;
		std::size_t BitsPerSample_;
		std::size_t BlockAlign_;
		std::size_t nChannels_;
		std::size_t AvgBytesPerSec_;
		std::size_t BufferLength_;
		std::size_t DataLength_;

		std::ofstream outputAudio_;

	public:
		encoder(void)
			: SamplePerSec_(0), BitsPerSample_(0), BlockAlign_(0),
			nChannels_(0), AvgBytesPerSec_(0), DataLength_(0)
		{

		}
		int open(const std::string &filename,
			std::size_t DataLength,
			std::size_t SamplePerSec = 48000,
			std::size_t BitsPerSample = 16)
		{
			// audio information
			SamplePerSec_ = SamplePerSec;
			BitsPerSample_ = BitsPerSample;
			DataLength_ = DataLength;
			nChannels_ = 1;
			BlockAlign_ = BitsPerSample_ / 8 * nChannels_;
			AvgBytesPerSec_ = SamplePerSec_ * BlockAlign_;
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

			fmt_samples_per_sec = SamplePerSec_;  // sampling frequency
			fmt_quantizationBitRate_per_sample = BitsPerSample_; // quantinzation bit

			data_chunk_ID[0] = 'd';
			data_chunk_ID[1] = 'a';
			data_chunk_ID[2] = 't';
			data_chunk_ID[3] = 'a';

			riff_chunk_size = 36 + DataLength_ * 2;
			fmt_channel = 1;
			fmt_bytes_per_sec = SamplePerSec_ * BitsPerSample_ / 8;
			fmt_block_size = BitsPerSample_ / 8;
			data_chunk_size = DataLength_ * 2;

			outputAudio_.write(riff_chunk_ID, 4);
			outputAudio_.write((char*)&riff_chunk_size, 4);
			outputAudio_.write(riff_form_type, 4);
			outputAudio_.write(fmt_chunk_ID, 4);
			outputAudio_.write((char*)&fmt_chunk_size, 4);
			outputAudio_.write((char*)&fmt_wave_format_type, 2);
			outputAudio_.write((char*)&fmt_channel, 2);	//
			outputAudio_.write((char*)&fmt_samples_per_sec, 4);
			outputAudio_.write((char*)&fmt_bytes_per_sec, 4);
			outputAudio_.write((char*)&fmt_block_size, 2);
			outputAudio_.write((char*)&fmt_quantizationBitRate_per_sample, 2);
			outputAudio_.write(data_chunk_ID, 4);
			outputAudio_.write((char*)&data_chunk_size, 4);

			return 0;
		}
		void at(std::size_t n, double value)
		{
			short data;
			double s;
			s = (value + 1.0) / 2.0 * 65536.0;
			if (s > 65536.0) s = 65536.0;
			if (s < 0.0) s = 0.0;
			data = static_cast<short>((s + 0.5) - 32768.0);
			outputAudio_.write(reinterpret_cast<char*>(&data), sizeof(short));
		}


	};

};

#endif

#endif // _AUDIO_HEADER_