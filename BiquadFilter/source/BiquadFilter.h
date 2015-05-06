//BiquadFilter.h
#pragma once

#include <vector>

namespace BiquadFilter
{
	//ディジタルフィルタ抽象クラス
	class Filter{
	public:
		virtual ~Filter();
		void filtering(std::vector<double> *x);
	protected:
		void alloc();

	protected:
		//フィルタ係数
		std::vector<double> a;
		std::vector<double> b;
	};

	//ローパスフィルタ
	class LPFilter : public Filter{
	public:
		LPFilter(double cutoff, double Q);
	private:
		double cutoff; //正規化されたカットオフ周波数(freq/samplerate)
		double Q;
	};

	//ハイパスフィルタ
	class HPFilter : public Filter{
	public:
		HPFilter(double cutoff, double Q);
	private:
		double cutoff; //正規化されたカットオフ周波数(freq/samplerate)
		double Q;
	};

	//バンドパスフィルタ
	class BPFilter : public Filter{
	public:
		BPFilter(double low_edge, double high_edge);
	private:
		double low_edge;
		double high_edge;
	};

	//ノッチフィルター
	class NTFilter : public Filter{
	public:
		NTFilter(double low_edge, double high_edge);
	private:
		double low_edge;
		double high_edge;
	};

	//ローシェルフフィルタ
	class LSFilter : public Filter{
	public:
		LSFilter(double cutoff, double Q, double gain);
	private:
		double cutoff;
		double Q;
		double gain;
	};

	//ハイシェルフフィルタ
	class HSFilter : public Filter{
	public:
		HSFilter(double cutoff, double Q, double gain);
	private:
		double cutoff;
		double Q;
		double gain;
	};

	//ピーキングフィルタ
	class PKFilter : public Filter{
	public:
		PKFilter(double low_edge, double high_edge, double gain);
	private:
		double low_edge;
		double high_edge;
		double gain;
	};

	//オールパスフィルタ
	class APFilter : public Filter{
	public:
		APFilter(double cutoff, double Q);
	private:
		double cutoff;
		double Q;
	};
}