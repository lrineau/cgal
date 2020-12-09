#undef QT_NO_KEYWORDS
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

//because dtk uses qt's foreach keyword, which we do not allow
#define foreach Q_FOREACH
#define signals Q_SIGNALS
#define slots Q_SLOTS

#include <dtkCore>
#include <dtkContinuousGeometry>
#include <dtkContinuousGeometryUtils>
#include <dtkBRep>
#include <dtkBRepReader>
#include <dtkNurbsSurface>

#include <QFileInfo>
#include <QDebug>

#include "Scene_cad_item.h"
#include "Messages_interface.h"
#include "Viewer.h"
#include "Scene.h"
#include "MainWindow.h"

class Cad_meshing_plugin :
    public QObject,
    public CGAL::Three::Polyhedron_demo_plugin_interface,
    public CGAL::Three::Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0")
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:
  //PLUGIN PART (for access to the scene and the message interface)
  bool applicable(QAction*) const
   {
     return false;
   }
   QList<QAction*> actions() const
   {
     return QList<QAction*>();
   }
   void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* sc, Messages_interface* mi) Q_DECL_OVERRIDE
   {
     this->messageInterface = mi;
     this->scene = sc;
     this->mw = mainWindow;
   }


  //IO PLUGIN PART
  QString name() const { return "Cad_meshing_plugin"; }
  QString nameFilters() const { return "CAD Files (*.3dm *.step *.stp *.iges *.igs)"; }
  bool canLoad(QFileInfo) const { return true; }
  QList<CGAL::Three::Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true) {
      auto ext = fileinfo.suffix().toLower();
      bool occt_file = (ext == "step") || (ext == "stp") || (ext == "igs") || (ext == "iges");
      if(!occt_file && ext != "3dm"){
        ok = false;
        return QList<CGAL::Three::Scene_item*>();
      }

      dtkLogger::instance().setLevel(dtk::LogLevel::Info);

      dtkContinuousGeometrySettings settings;
      settings.beginGroup("continuous-geometry");

      dtkContinuousGeometry::setVerboseLoading(true);
      QString plugins_path = settings.value("plugins").toString();
      if(plugins_path.isEmpty())
      {
        plugins_path = QFileDialog::getExistingDirectory(mw, "Path to install/plugins/dtkContinuousPlugins");
      }
      dtkContinuousGeometry::initialize(plugins_path);
      settings.endGroup();
      dtkBRepReader *brep_reader =
          dtkContinuousGeometry::bRepReader::pluginFactory().create(
              occt_file ? "OpenCASCADEBRepReader" : "openNURBSBRepReader");
      if (brep_reader == NULL) {
        messageInterface->message_error("ERROR in initialization");
        ok = false;
        return QList<CGAL::Three::Scene_item*>();
      }
      brep_reader->setInputBRepFilePath(fileinfo.absoluteFilePath());
      brep_reader->run();
      dtkBRep* brep = brep_reader->outputDtkBRep();
      Scene_cad_item* item = new Scene_cad_item(brep, scene);
      item->setName(fileinfo.baseName());
      item->setFlatMode();
      QList<int> children_ids;
      for( CGAL::Three::Scene_interface::Item_id id : item->getChildren())
      {
        CGAL::Three::Scene_item* child =scene->item(id);
        children_ids.append(scene->item_id(child));
      }
      static_cast<Scene*>(scene)->setSelectedItemsList(children_ids);
      static_cast<MainWindow*>(this->mw)->colorItems();
      scene->setSelectedItem(scene->item_id(item));
      if(add_to_scene)
        CGAL::Three::Three::scene()->addItem(item);
      ok = true;
      return QList<CGAL::Three::Scene_item*>()<<item;
  }

  bool canSave(const CGAL::Three::Scene_item*) {
      return false;
  }

  bool save(QFileInfo, QList<CGAL::Three::Scene_item*>& ){
    return false;
  }


  Messages_interface* messageInterface;
  //The reference to the scene
  CGAL::Three::Scene_interface* scene;
  //The reference to the main window
  QMainWindow* mw;
};
#include "Cad_meshing_plugin.moc"
