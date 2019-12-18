/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// VertReg Logic includes
#include <vtkSlicerVertRegLogic.h>

// VertReg includes
#include "qSlicerVertRegModule.h"
#include "qSlicerVertRegModuleWidget.h"

//-----------------------------------------------------------------------------
#if (QT_VERSION < QT_VERSION_CHECK(5, 0, 0))
#include <QtPlugin>
Q_EXPORT_PLUGIN2(qSlicerVertRegModule, qSlicerVertRegModule);
#endif

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerVertRegModulePrivate
{
public:
  qSlicerVertRegModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerVertRegModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerVertRegModulePrivate::qSlicerVertRegModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerVertRegModule methods

//-----------------------------------------------------------------------------
qSlicerVertRegModule::qSlicerVertRegModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerVertRegModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerVertRegModule::~qSlicerVertRegModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerVertRegModule::helpText() const
{
  return "This is a loadable module that can be bundled in an extension";
}

//-----------------------------------------------------------------------------
QString qSlicerVertRegModule::acknowledgementText() const
{
  return "This work was partially funded by NIH grant NXNNXXNNNNNN-NNXN";
}

//-----------------------------------------------------------------------------
QStringList qSlicerVertRegModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString("John Doe (AnyWare Corp.)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerVertRegModule::icon() const
{
  return QIcon(":/Icons/VertReg.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerVertRegModule::categories() const
{
  return QStringList() << "Examples";
}

//-----------------------------------------------------------------------------
QStringList qSlicerVertRegModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerVertRegModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerVertRegModule
::createWidgetRepresentation()
{
  return new qSlicerVertRegModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerVertRegModule::createLogic()
{
  return vtkSlicerVertRegLogic::New();
}
