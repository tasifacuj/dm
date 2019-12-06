#include "MainWindow.hpp"
#include "kd_tree/PointCloud.hpp"
#include "kd_tree/kd-tree.hpp"

#include <QDebug>

#include <cstdlib>

namespace dm {


template <typename T>
void generateRandomPointCloud(KdTree::PointCloud<T> &point, const size_t N, const T max_range = 10)
{
	// Generating Random Point Cloud
	point.pts.resize(N);
	for (size_t i = 0; i < N; i++)
	{
		point.pts[i].x = max_range * (rand() % 1000) / T(1000);
		point.pts[i].y = max_range * (rand() % 1000) / T(1000);
		point.pts[i].z = max_range * (rand() % 1000) / T(1000);
	}
}

template <typename num_t>
void kdtree_demo(const size_t N)
{
	KdTree::PointCloud<num_t> cloud;

	// Generate points:
	generateRandomPointCloud(cloud, N);

	// construct a kd-tree index:
	typedef KdTree::KDTreeSingleIndexAdaptor<
		KdTree::L2_Simple_Adaptor<num_t, KdTree::PointCloud<num_t> >,
		KdTree::PointCloud<num_t>,
		3 /* dim */
	> my_kd_tree_t;

	my_kd_tree_t   index(3 /*dim*/, cloud, KdTree::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	index.buildIndex();
	const num_t query_pt[3] = { 0.5, 0.5, 0.5 };
	// ----------------------------------------------------------------
	// radiusSearch(): Perform a search for the points within search_radius
	// ----------------------------------------------------------------
	{
		const num_t search_radius = static_cast<num_t>(0.1);
		std::vector<std::pair<size_t, num_t> >   ret_matches;

		KdTree::SearchParams params;
		//params.sorted = false;

		const size_t nMatches = index.radiusSearch(&query_pt[0], search_radius, ret_matches, params);

		qDebug() << "radiusSearch(): radius=" << search_radius << " -> " << nMatches << " matches";
		
		for (size_t i = 0; i < nMatches; i++)
			qDebug() << "idx[" << i << "]=" << ret_matches[i].first << " dist[" << i << "]=" << ret_matches[i].second;
	}
}

MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	QObject::connect(ui.btn_, &QPushButton::clicked, this, &MainWindow::onClickMeClicked);
	srand(static_cast<unsigned int>(time(nullptr)));
}


void MainWindow::onClickMeClicked() {

	//kdtree_demo<float>(20);
	kdtree_demo<double>(100000);
	qDebug() << __func__;
}

}