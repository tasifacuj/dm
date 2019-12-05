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
}

MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	QObject::connect(ui.btn_, &QPushButton::clicked, this, &MainWindow::onClickMeClicked);
}


void MainWindow::onClickMeClicked() {
	kdtree_demo<float>(4);
	qDebug() << __func__;
}

}