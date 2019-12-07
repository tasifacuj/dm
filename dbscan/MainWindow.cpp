#include "MainWindow.hpp"
//#include "kd_tree/PointCloud.hpp"
#include "kd_tree/kd-tree.hpp"
#include "dbscan/PointCloud2D.hpp"

#include <QDebug>

#include <cstdlib>

namespace dm {


template <typename T>
void generateRandomPointCloud(dbscan::PointCloud2D<T> &point, const size_t N, const T max_range = 10)
{
	// Generating Random Point Cloud
	point.pts.resize(N);
	
	for (size_t i = 0; i < N; i++){
		point.pts[i].x = max_range * (rand() % 1000) / T(1000);
		point.pts[i].y = max_range * (rand() % 1000) / T(1000);
		//point.pts[i].z = max_range * (rand() % 1000) / T(1000);
	}
}

template <typename num_t>
void kdtree_demo(const size_t N)
{
	dbscan::PointCloud2D<num_t> cloud;

	// Generate points:
	generateRandomPointCloud(cloud, N);

	// construct a kd-tree index:
	typedef KdTree::KDTreeSingleIndexAdaptor<
		KdTree::L2_Simple_Adaptor<num_t, dbscan::PointCloud2D<num_t> >,
		dbscan::PointCloud2D<num_t>,
		2 /* dim */
	> my_kd_tree_t;

	my_kd_tree_t   index(2 /*dim*/, cloud, KdTree::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	index.buildIndex();
	const num_t query_pt[2] = { 0.5, 0.5 };
	// ----------------------------------------------------------------
	// radiusSearch(): Perform a search for the points within search_radius
	// ----------------------------------------------------------------
	{
		const num_t search_radius = static_cast<num_t>(0.01);
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