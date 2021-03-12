#undef QT_NO_KEYWORDS
#define foreach Q_FOREACH
#define signals Q_SIGNALS
#define slots Q_SLOTS
#include <set>
#include <dtkBRep>
#include <dtkNurbsSurface>
#include <dtkTopoTrim>
#include <dtkNurbsCurve>
#include <dtkTrimLoop>
#include <dtkTrim>
#include <dtkContinuousGeometryUtils>

#include <CGAL/Three/Edge_container.h>
#include <CGAL/Three/Three.h>

#include "Scene_cad_item.h"
#include "Scene_nurbs_item.h"
#include <QDebug>
#include <QMenu>


typedef CGAL::Three::Edge_container Ec;
typedef CGAL::Three::Viewer_interface Vi;

struct Scene_cad_item_priv{
  Scene_cad_item_priv(dtkBRep* brep, CGAL::Three::Scene_interface* scene, Scene_cad_item* parent)
    :m_brep(brep), item(parent), intersections_shown(false)
  {
    std::size_t i = 0;
    const std::vector < dtkNurbsSurface* > nurbs_surfaces = m_brep->nurbsSurfaces();

    for (auto it = nurbs_surfaces.begin(); it != nurbs_surfaces.end(); ++it) {
      Scene_nurbs_item* nurbs_item =  new Scene_nurbs_item(*(*it));
      nurbs_item->setFlatMode();
      nurbs_item->setName(QString("NURBS Surface #%1").arg(i));
      scene->addItem(nurbs_item);
      item->addChild(nurbs_item);
      nurbs_item->moveToGroup(item);
      item->connect(item, SIGNAL(highlighted(const dtkTopoTrim *)), nurbs_item, SLOT(highlight(const dtkTopoTrim *)));
      ++i;
    }

    std::set< const dtkTopoTrim *> topo_trims;

    intersection.resize(0);
    intersection_colors.resize(0);

    std::size_t index = 0;
    for(auto topo_trim : brep->topoTrims()) {
        intersection_colors_indices.insert(topo_trim, index);
        color_indices.insert(index);
        dtkNurbsCurve* nurb = topo_trim->m_nurbs_curve_3d;
        if (nurb != nullptr) {
            std::size_t length = nurb->nbCps() + nurb->degree() - 1;
            std::vector<double> knots(length);
            nurb->knots(knots.data());
            dtkContinuousGeometryPrimitives::Point_3 p(0,0,0);
            nurb->evaluatePoint(knots[0], p.data());
            for(int j=0; j<3; ++j)
                intersection.push_back(p[j]);
            ++index;
            for(float f = knots[0]+1/100.0*(knots[length-1]-knots[0]);
                f<knots[length-1] - 1/100.0*(knots[length-1]-knots[0]);
                f+=1/100.0*(knots[length-1]-knots[0]))
                {
                    nurb->evaluatePoint(f, p.data());
                    for(int j=0; j<3; ++j) {
                        intersection.push_back(p[j]);
                    }
                    intersection_colors.push_back((float)0);
                    intersection_colors.push_back((float)0);
                    intersection_colors.push_back((float)0);
                    ++index;
                    for(int j=0; j<3; ++j) {
                        intersection.push_back(p[j]);
                    }
                    intersection_colors.push_back((float)0);
                    intersection_colors.push_back((float)0);
                    intersection_colors.push_back((float)0);
                    ++index;
                }

            nurb->evaluatePoint(knots[length-1], p.data());
            for(int j=0; j<3; ++j)
                intersection.push_back(p[j]);
            ++index;
            intersection_colors.push_back((float)0);
            intersection_colors.push_back((float)0);
            intersection_colors.push_back((float)0);
        }
        }
  }


    dtkBRep* m_brep;

    Scene_cad_item* item;

    mutable QOpenGLShaderProgram* m_program;

    mutable bool intersections_shown;
    mutable std::vector<float> intersection;
    mutable QHash<const dtkTopoTrim *, std::size_t> intersection_colors_indices; //starts from 1
    mutable std::set< std::size_t > color_indices;

    mutable std::size_t current_index;
    mutable std::size_t current_next_index;
    mutable std::vector<float> intersection_colors;

    enum Vao
        {
            INTERSECTION=0,
            NumberOfVaos
        };

    enum Buffer
        {
            B_INTERSECTION=0,
            B_INTERSECTION_COLORS=1,
            NumberOfBuffers
        };
};


Scene_cad_item::Scene_cad_item(dtkBRep* brep, CGAL::Three::Scene_interface* scene)
  :CGAL::Three::Scene_group_item("unnamed")
{
  setEdgeContainer(0,
                   new Ec( Vi::PROGRAM_NO_SELECTION
                           ,false));
  d = new Scene_cad_item_priv(brep, scene, this);
  d->current_index = 0;
  d->current_next_index = 0;
}

void Scene_cad_item::computeElements()const
{
  getEdgeContainer(0)->allocate(Ec::Vertices, d->intersection.data(),
                                          static_cast<int>(d->intersection.size() *
                                                           sizeof(float)));
  getEdgeContainer(0)->allocate(Ec::Colors, d->intersection_colors.data(),
                                          static_cast<int>(d->intersection_colors.size() *
                                                           sizeof(float)));

}


void Scene_cad_item::draw(CGAL::Three::Viewer_interface* viewer) const
{
  CGAL::Three::Scene_group_item::draw(viewer);
  if(!isInit(viewer))
    initGL(viewer);
  if ( getBuffersFilled() &&
       ! getBuffersInit(viewer))
  {
    initializeBuffers(viewer);
    setBuffersInit(viewer, true);
  }
  if(!getBuffersFilled())
  {
    initializeBuffers(viewer);
    setBuffersFilled(true);
    setBuffersInit(viewer, true);
  }

  if(d->intersections_shown)
  {
    getEdgeContainer(0)->setColor(QColor(Qt::red));
    getEdgeContainer(0)->draw(viewer, false);
  }
}

void Scene_cad_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const
{
  CGAL::Three::Scene_group_item::drawEdges(viewer);
}
Scene_cad_item::~Scene_cad_item()
{
  if(d)
    delete d;
}

Scene_cad_item::Bbox Scene_cad_item::bbox() const
{
  Bbox box = Bbox(0,0,0,0,0,0);
  for(CGAL::Three::Scene_interface::Item_id id: getChildren())
  {
    Scene_item* item = scene->item(id);
    box += item->bbox();
  }
  return box;
}

void Scene_cad_item::show_trimmed(bool b)
{
  for(CGAL::Three::Scene_interface::Item_id id: getChildren())
  {
    Scene_item* item = scene->item(id);
    Scene_nurbs_item* nurbs = qobject_cast<Scene_nurbs_item*>(item);
    if(!nurbs)
      continue;
    QAction* actionShowTrimmed = nurbs->contextMenu()->findChild<QAction*>("actionShowTrimmed");
    if(!actionShowTrimmed){
      continue;
    }
    actionShowTrimmed->toggled(!b);
    //nurbs->show_trimmed(b);
  }
}

void Scene_cad_item::show_control_points(bool b)
{
  for(CGAL::Three::Scene_interface::Item_id id: getChildren())
  {
    Scene_item* item = scene->item(id);
    Scene_nurbs_item* nurbs = qobject_cast<Scene_nurbs_item*>(item);
    if(!nurbs)
      continue;
    QAction* actionShowCPs= nurbs->contextMenu()->findChild<QAction*>("actionShowCPs");
    if(!actionShowCPs){
      continue;
    }
    actionShowCPs->toggled(!b);
  }
}

void Scene_cad_item::show_bezier_surfaces(bool b)
{
  for(CGAL::Three::Scene_interface::Item_id id: getChildren())
  {
    Scene_item* item = scene->item(id);
    Scene_nurbs_item* nurbs = qobject_cast<Scene_nurbs_item*>(item);
    if(!nurbs)
      continue;
    QAction* actionShowBeziers = nurbs->contextMenu()->findChild<QAction*>("actionShowBeziers");
    if(!actionShowBeziers){
      continue;
    }
    actionShowBeziers->toggled(!b);
    nurbs->setVisible(!b);
  }
}

void Scene_cad_item::show_intersections(bool b)
{
  d->intersections_shown = b;
  itemChanged();
}

QMenu* Scene_cad_item::contextMenu()
{
  const char* prop_name = "Menu modified by Scene_cad_item.";

  QMenu* menu = Scene_item::contextMenu();

  // Use dynamic properties:
  // http://doc.qt.io/qt-5/qobject.html#property
  bool menuChanged = menu->property(prop_name).toBool();

  if (!menuChanged) {
    QAction* actionShowTrimmed=
        menu->addAction(tr("Show Trimmed"));
    actionShowTrimmed->setCheckable(true);
    actionShowTrimmed->setChecked(false);
    actionShowTrimmed->setObjectName("actionShowTrimmed");
    actionShowTrimmed->setProperty("is_groupable", true);
    connect(actionShowTrimmed, SIGNAL(toggled(bool)),
            this, SLOT(show_trimmed(bool)));

    QAction* actionShowCPs=
        menu->addAction(tr("Show Control Points"));
    actionShowCPs->setCheckable(true);
    actionShowCPs->setChecked(false);
    actionShowCPs->setProperty("is_groupable", true);
    actionShowCPs->setObjectName("actionShowCPs");
    connect(actionShowCPs, SIGNAL(toggled(bool)),
            this, SLOT(show_control_points(bool)));

    QAction* actionShowBezierSurfaces=
        menu->addAction(tr("Show Bezier Surfaces"));
    actionShowBezierSurfaces->setCheckable(true);
    actionShowBezierSurfaces->setChecked(false);
    actionShowBezierSurfaces->setProperty("is_groupable", true);
    actionShowBezierSurfaces->setObjectName("actionShowBezierSurfaces");
    connect(actionShowBezierSurfaces, SIGNAL(toggled(bool)),
            this, SLOT(show_bezier_surfaces(bool)));

    QAction* actionShowIntersections=
        menu->addAction(tr("Show Intersections"));
    actionShowIntersections->setCheckable(true);
    actionShowIntersections->setChecked(false);
    actionShowIntersections->setObjectName("actionShowIntersections");
    actionShowIntersections->setProperty("is_groupable", true);
    connect(actionShowIntersections, SIGNAL(toggled(bool)),
            this, SLOT(show_intersections(bool)));

    menu->setProperty(prop_name, true);
  }
  return menu;
}

dtkBRep* Scene_cad_item::brep()
{
  return d->m_brep;
}

const dtkBRep* Scene_cad_item::brep() const
{
  return d->m_brep;
}

QString Scene_cad_item::toolTip() const
{
  QString str =
         QObject::tr("<p>BRep <b>%1</b>")
            .arg(this->name());
  return str;
}

void Scene_cad_item::clearHighlight(void)
{
    if(d->current_next_index !=0) {
        std::vector<float> old_intersection_colors((d->current_next_index - d->current_index) * 3, 0.);

        getEdgeContainer(0)->getVbo(Ec::Colors)->bind();
        getEdgeContainer(0)->getVbo(Ec::Colors)->vbo.write(d->current_index * 3 * sizeof(float), old_intersection_colors.data(), 3 * (d->current_next_index - d->current_index) * sizeof(float));
        getEdgeContainer(0)->getVbo(Ec::Colors)->release();
    }

    d->current_index = 0;
    d->current_next_index = 0;

    itemChanged();
}

void Scene_cad_item::highlight(const dtkTopoTrim *topo_trim)
{
    // ///////////////////////////////////////////////////////////////////
    // Clear the buffers
    // ///////////////////////////////////////////////////////////////////
    if(d->current_next_index !=0) {
        std::vector<float> old_intersection_colors((d->current_next_index - d->current_index) * 3, 0.);
         getEdgeContainer(0)->getVbo(Ec::Colors)->bind();
         getEdgeContainer(0)->getVbo(Ec::Colors)->vbo.write(d->current_index * 3 * sizeof(float), old_intersection_colors.data(), 3 * (d->current_next_index - d->current_index) * sizeof(float));
         getEdgeContainer(0)->getVbo(Ec::Colors)->release();

    }

    // ///////////////////////////////////////////////////////////////////
    // Recovers the index associated to the trim
    // ///////////////////////////////////////////////////////////////////
    std::size_t index = 0;
    auto it = d->intersection_colors_indices.find(topo_trim);
    if(it != d->intersection_colors_indices.end()) {
        index = it.value();
    } else {
        dtkWarn() << "The pointer to the dtkTopoTrim is not listed in the possible trims to display.";
        return;
    }

    // ///////////////////////////////////////////////////////////////////
    // Finds next index to count the number of spheres to change
    // ///////////////////////////////////////////////////////////////////
    auto next = d->color_indices.upper_bound(index);
    std::size_t next_index = 0;
    if(next != d->color_indices.end()) {
        next_index = *next;
    } else {
        next_index = d->intersection_colors.size() / 3;
    }

    std::vector<float> new_intersection_colors((next_index - index) * 3);
    for(std::size_t i = 0; i < next_index - index; ++i) {
        new_intersection_colors[3 * i] =     (float)255;
        new_intersection_colors[3 * i + 1] = (float)255;
        new_intersection_colors[3 * i + 2] = (float)0;
    }

    getEdgeContainer(0)->getVbo(Ec::Colors)->bind();
    getEdgeContainer(0)->getVbo(Ec::Colors)->vbo.write(index * 3 * sizeof(float), new_intersection_colors.data(), 3 * (next_index - index) * sizeof(float));
    getEdgeContainer(0)->getVbo(Ec::Colors)->release();


    d->current_index = index;
    d->current_next_index = next_index;
    itemChanged();
    Q_EMIT highlighted(topo_trim);
    Q_EMIT updated();
}

void Scene_cad_item::initializeBuffers(Viewer_interface *viewer) const
{

  getEdgeContainer(0)->initializeBuffers(viewer);
  getEdgeContainer(0)->setFlatDataSize(d->intersection.size()/3);

  d->intersection.clear();
  d->intersection.shrink_to_fit();
}
