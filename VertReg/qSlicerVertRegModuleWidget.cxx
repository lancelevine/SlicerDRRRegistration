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

// Qt includes
#include <QDebug>

// SlicerQt includes
#include "qSlicerVertRegModuleWidget.h"
#include "ui_qSlicerVertRegModuleWidget.h"

#include "qSlicerCoreIOManager.h"

#include "vtkSlicerApplicationLogic.h"
#include "vtkMRMLSegmentationNode.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCommand.h"
#include "itkSimpleFilterWatcher.h"

#include "itkMeanProjectionImageFilter.h"
#include "itkProjectionImageFilter.h"


#include "itkMRMLIDImageIO.h"
#include "itkMRMLIDImageIOFactory.h"
#include "itkVTKImageToImageFilter.h"
#include "itkImageToVTKImageFilter.h"

#include "itkVTKImageIO.h"

#include "vtkMatrix4x4.h"
#include "itkFlipImageFilter.h"

#include "vtkRenderer.h"
#include "vtkImageActor.h"
#include "vtkJPEGReader.h"

#include "itkResampleImageFilter.h"
#include "itkCenteredEuler3DTransform.h"
#include "itkRayCastInterpolateImageFunction.h"
#include "itkRescaleIntensityImageFilter.h"
#include "vtkImageReslice.h"

#include "vtkImageWriter.h"
#include "vtkNrrdReader.h"
#include "itkExtractImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkPasteImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkNrrdImageIO.h"
#include "itkVTKImageIO.h"

#include "vtkImageCast.h"
#include "vtkPNGWriter.h"
#include "vtkImageShiftScale.h"

#include "itkCenteredAffineTransform.h"

#include "itkEuler2DTransform.h"
#include "itkExhaustiveOptimizerv4.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkCenteredTransformInitializer.h"
#include "itkImageRegistrationMethodv4.h"

#include "itkVersorRigid3DTransformOptimizer.h"
#include "itkVersorRigid3DTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkCommand.h"

#include "itkRescaleIntensityImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkImageRegistrationMethodv4.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"
#include "itkCenteredTransformInitializer.h"
#include "itkSimilarity2DTransform.h"

#include <vtkTransform.h>
#include "itkSubtractImageFilter.h"

#include <chrono>


class vtkMRMLNode;

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerVertRegModuleWidgetPrivate: public Ui_qSlicerVertRegModuleWidget
{
public:
  qSlicerVertRegModuleWidgetPrivate();

  vtkWeakPointer<vtkMRMLVolumeNode> InputVolumeNode;
  vtkWeakPointer<vtkMRMLVolumeNode> OutputVolumeNode;
};

//-----------------------------------------------------------------------------
// qSlicerVertRegModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerVertRegModuleWidgetPrivate::qSlicerVertRegModuleWidgetPrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerVertRegModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerVertRegModuleWidget::qSlicerVertRegModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerVertRegModuleWidgetPrivate )
{
}

//-----------------------------------------------------------------------------
qSlicerVertRegModuleWidget::~qSlicerVertRegModuleWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerVertRegModuleWidget::setup()
{
  Q_D(qSlicerVertRegModuleWidget);
  d->setupUi(this);
  this->Superclass::setup();

  debug = true;

  connect(d->InputVolumeComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode*)),
	  this, SLOT(setInputVolume(vtkMRMLNode*)));

  connect(d->ProjectionComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode*)),
	  this, SLOT(setOutputVolume(vtkMRMLNode*)));

  connect(d->ProjectionButton0, SIGNAL(clicked()),
	  this, SLOT(projectVolumeDRR0()));

  connect(d->ProjectionButton1, SIGNAL(clicked()),
	  this, SLOT(projectVolumeDRR5()));

  connect(d->ProjectionButton11, SIGNAL(clicked()),
	  this, SLOT(projectVolumeDRR15()));

  connect(d->ProjectionButton12, SIGNAL(clicked()),
	  this, SLOT(projectVolumeDRR25()));

  connect(d->ProjectionButton2, SIGNAL(clicked()),
	  this, SLOT(rotateAndProjectVolume0()));

  connect(d->ProjectionButton3, SIGNAL(clicked()),
	  this, SLOT(rotateAndProjectVolume5()));

  connect(d->ProjectionButton4, SIGNAL(clicked()),
	  this, SLOT(rotateAndProjectVolume15()));

  connect(d->ProjectionButton5, SIGNAL(clicked()),
	  this, SLOT(rotateAndProjectVolume25()));

  connect(d->RegisterButton, SIGNAL(clicked()),
	  this, SLOT(registerVolume2dSimilarity()));

  connect(d->RegisterButton2, SIGNAL(clicked()),
	  this, SLOT(registerVolume2dEuler()));

  connect(d->RegisterButton3, SIGNAL(clicked()),
	  this, SLOT(automaticRegistration()));

  connect(d->CropButton, SIGNAL(clicked()),
	  this, SLOT(rotateVolume()));

}

/*template <typename ITK_Exporter, typename VTK_Importer>
void ConnectPipelines(ITK_Exporter exporter, VTK_Importer* importer)
{
	importer->SetUpdateInformationCallback(exporter->GetUpdateInformationCallback());
	importer->SetPipelineModifiedCallback(exporter->GetPipelineModifiedCallback());
	importer->SetWholeExtentCallback(exporter->GetWholeExtentCallback());
	importer->SetSpacingCallback(exporter->GetSpacingCallback());
	importer->SetOriginCallback(exporter->GetOriginCallback());
	importer->SetScalarTypeCallback(exporter->GetScalarTypeCallback());
	importer->SetNumberOfComponentsCallback(exporter->GetNumberOfComponentsCallback());
	importer->SetPropagateUpdateExtentCallback(exporter->GetPropagateUpdateExtentCallback());
	importer->SetUpdateDataCallback(exporter->GetUpdateDataCallback());
	importer->SetDataExtentCallback(exporter->GetDataExtentCallback());
	importer->SetBufferPointerCallback(exporter->GetBufferPointerCallback());
	importer->SetCallbackUserData(exporter->GetCallbackUserData());
}*/

//-----------------------------------------------------------------------------
void qSlicerVertRegModuleWidget::setInputVolume(vtkMRMLNode* volumeNode)
{
	Q_D(qSlicerVertRegModuleWidget);
	qDebug() << "function::setInputVolume";

	//qvtkReconnect(d->InputVolumeNode, volumeNode, vtkCommand::ModifiedEvent, this, SLOT(updateVolumeInfo()));
	d->InputVolumeNode = volumeNode;
	//qDebug() << volumeNode->GetID();
}

//-----------------------------------------------------------------------------
void qSlicerVertRegModuleWidget::setOutputVolume(vtkMRMLNode* volumeNode)
{
	Q_D(qSlicerVertRegModuleWidget);
	qDebug() << "function::setOutputVolume";

	//qvtkReconnect(d->InputVolumeNode, volumeNode, vtkCommand::ModifiedEvent, this, SLOT(updateVolumeInfo()));
	d->OutputVolumeNode = volumeNode;
	//qDebug() << volumeNode->GetID();
}

void qSlicerVertRegModuleWidget::projectVolumeFilter()
{
	Q_D(qSlicerVertRegModuleWidget);

	using clock = std::chrono::system_clock;
	using ms = std::chrono::milliseconds;

	const auto before = clock::now();

	qDebug() << "function::projectVolumeFilter()";

	const int dim = 3;

	typedef short PType;
	typedef itk::Image< PType, dim > InputImageType;
	//typedef itk::Image< PType, dim > OutputImageType2;
	typedef itk::Image< PType, dim - 1 > OutputImageType;


	vtkSmartPointer<vtkImageData> inputVTK = vtkSmartPointer<vtkImageData>::New();
	inputVTK = d->InputVolumeNode->GetImageData();

	//vtkNew<vtkMatrix4x4> matr;
	//d->InputVolumeNode->GetIJKToRASMatrix(matr);
	//d->OutputVolumeNode->SetIJKToRASMatrix(matr);

	qDebug() << "function::projectVolumeFilter1()";

	/*vtkSmartPointer<vtkImageCast> castFilter =
	vtkSmartPointer<vtkImageCast>::New();
	castFilter->SetOutputScalarTypeToUnsignedChar();
	castFilter->SetInputData(inputVTK);
	castFilter->Update();*/

	qDebug() << "function::projectVolumeFilter12()";


	typedef itk::VTKImageToImageFilter< InputImageType > VTKtoITKFilterType;
	VTKtoITKFilterType::Pointer vtkToItkFilter = VTKtoITKFilterType::New();
	vtkToItkFilter->SetInput(inputVTK);
	vtkToItkFilter->Update();
	InputImageType::Pointer image = vtkToItkFilter->GetOutput();
	
	qDebug() << "function::projectVolumeFilter2()";

	using FilterType = itk::ThresholdImageFilter< InputImageType >;
	FilterType::Pointer thresholdfilter = FilterType::New();
	thresholdfilter->SetInput(image);
	thresholdfilter->ThresholdOutside(95, 500);
	thresholdfilter->SetOutsideValue(0);
	thresholdfilter->Update();


	typedef itk::MeanProjectionImageFilter< InputImageType, OutputImageType > ProjectionFilterType;
	ProjectionFilterType::Pointer projfilter = ProjectionFilterType::New();
	projfilter->SetInput(thresholdfilter->GetOutput());
	projfilter->SetProjectionDimension(1);
	projfilter->Update();

	//itk::SimpleFilterWatcher watcher(projfilter, "filter");
	qDebug() << "function::projectVolumeFilter3()";


	OutputImageType::Pointer outputITK = projfilter->GetOutput();

	qDebug() << "function::projectVolumeFilter4()";


	if (debug) {
		/*typedef itk::Image< unsigned char, dim - 1 > WriterImageType;
		using CastFilterType = itk::CastImageFilter<OutputImageType, WriterImageType >;
		CastFilterType::Pointer castFilteritk = CastFilterType::New();
		castFilteritk->SetInput(outputITK);
		castFilteritk->Update();

		typedef itk::ImageFileWriter< WriterImageType > WriterType;*/
		typedef itk::ImageFileWriter< OutputImageType > WriterType; 
		WriterType::Pointer writer = WriterType::New();
		//writer->SetInput(castFilteritk->GetOutput());
		writer->SetInput(outputITK);
		writer->SetFileName("output_moving.nrrd");
		writer->Update();
	}

	qDebug() << "function::projectVolumeFilter5()";

	using ReaderType = itk::ImageFileReader< OutputImageType >;
	ReaderType::Pointer reader2 = ReaderType::New();
	reader2->SetFileName("output_moving.nrrd");
	reader2->Update();

	/*using ITKtoVTKFilterType2 = itk::ImageToVTKImageFilter< OutputImageType >;
	ITKtoVTKFilterType2::Pointer itktovtkfilter2 = ITKtoVTKFilterType2::New();
	itktovtkfilter2->SetInput(reader2->GetOutput());
	itktovtkfilter2->Update();*/
	vtkSmartPointer<vtkImageData> outputVTK = vtkSmartPointer<vtkImageData>::New();
	//outputVTK = itktovtkfilter2->GetOutput();

	//outputVTK->Print(std::cout);
	//vtkImageData * outputVTK = itktovtkfilter->GetOutput();

	qDebug() << "function::projectVolumeFilter6()";

	if (debug) {

		vtkSmartPointer<vtkNrrdReader> writer = vtkSmartPointer<vtkNrrdReader>::New();
		writer->SetInputData(outputVTK);
		writer->SetFileName("output_moving.nrrd");
		writer->Update();
		outputVTK = writer->GetOutput();
	}


	/*vtkSmartPointer<vtkImageCast> castFilter2 =
		vtkSmartPointer<vtkImageCast>::New();
	castFilter2->SetOutputScalarTypeToShort();
	castFilter2->SetInputData(outputVTK);
	castFilter2->Update();*/

	//d->OutputVolumeNode->SetAndObserveImageData(castFilter->GetOutput());
	d->OutputVolumeNode->SetAndObserveImageData(outputVTK);

	const auto duration = std::chrono::duration_cast<ms>(clock::now() - before);

	std::cout << "It took " << duration.count() / 1000.0 << "ms" << std::endl;
}

//3d
/*class CommandIterationUpdate : public itk::Command
{
public:
	typedef  CommandIterationUpdate   Self;
	typedef  itk::Command             Superclass;
	typedef itk::SmartPointer<Self>   Pointer;
	itkNewMacro(Self);
protected:
	CommandIterationUpdate() {};
public:
	typedef itk::VersorRigid3DTransformOptimizer OptimizerType;
	typedef   const OptimizerType *              OptimizerPointer;
	void Execute(itk::Object *caller, const itk::EventObject & event)
	{
		Execute((const itk::Object *)caller, event);
	}
	void Execute(const itk::Object * object, const itk::EventObject & event)
	{
		OptimizerPointer optimizer = static_cast< OptimizerPointer >(object);
		if (!itk::IterationEvent().CheckEvent(&event))
		{
			return;
		}
		std::cout << optimizer->GetCurrentIteration() << "   ";
		std::cout << optimizer->GetValue() << "   ";
		std::cout << optimizer->GetCurrentPosition() << std::endl;
	}
};*/

//similarity
/*class CommandIterationUpdate : public itk::Command
{
public:
	using Self = CommandIterationUpdate;
	using Superclass = itk::Command;
	using Pointer = itk::SmartPointer<Self>;
	itkNewMacro(Self);
protected:
	CommandIterationUpdate() = default;
public:
	using OptimizerType = itk::RegularStepGradientDescentOptimizerv4<double>;
	using OptimizerPointer = const OptimizerType *;
	void Execute(itk::Object *caller, const itk::EventObject & event) override
	{
		Execute((const itk::Object *)caller, event);
	}
	void Execute(const itk::Object * object, const itk::EventObject & event) override
	{
		auto optimizer = static_cast< OptimizerPointer >(object);
		if (!itk::IterationEvent().CheckEvent(&event))
		{
			return;
		}
		std::cout << optimizer->GetCurrentIteration() << "   ";
		std::cout << optimizer->GetValue() << "   ";
		std::cout << optimizer->GetCurrentPosition() << std::endl;
	}
};*/

//euler
class CommandIterationUpdate : public itk::Command
{
public:
	using Self = CommandIterationUpdate;
	using Superclass = itk::Command;
	using Pointer = itk::SmartPointer<Self>;
	itkNewMacro(Self);
protected:
	CommandIterationUpdate() = default;
public:
	using OptimizerType = itk::RegularStepGradientDescentOptimizerv4<double>;
	using OptimizerPointer = const OptimizerType *;
	void Execute(itk::Object *caller, const itk::EventObject & event) override
	{
		Execute((const itk::Object *)caller, event);
	}
	void Execute(const itk::Object * object, const itk::EventObject & event) override
	{
		auto optimizer = static_cast< OptimizerPointer >(object);
		if (!itk::IterationEvent().CheckEvent(&event))
		{
			return;
		}
		std::cout << optimizer->GetCurrentIteration() << "   ";
		std::cout << optimizer->GetValue() << "   ";
		std::cout << optimizer->GetCurrentPosition() << std::endl;
	}
};

void qSlicerVertRegModuleWidget::registerVolume3d()
{
	Q_D(qSlicerVertRegModuleWidget);

	using clock = std::chrono::system_clock;
	using ms = std::chrono::milliseconds;

	const auto before = clock::now();

	qDebug() << "function::registerVolume()";

	const unsigned int                          Dimension = 3;
	typedef  float                              PixelType;
	typedef itk::Image< PixelType, Dimension >  FixedImageType;
	typedef itk::Image< PixelType, Dimension >  MovingImageType;
	//  Software Guide : BeginLatex
	//
	//  The Transform class is instantiated using the code below. The only
	//  template parameter to this class is the representation type of the
	//  space coordinates.
	//
	//  \index{itk::Versor\-Rigid3D\-Transform!Instantiation}
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	typedef itk::VersorRigid3DTransform< double > TransformType;
	// Software Guide : EndCodeSnippet
	typedef itk::VersorRigid3DTransformOptimizer                                  OptimizerType;
	typedef itk::MeanSquaresImageToImageMetric< FixedImageType, MovingImageType > MetricType;
	typedef itk::LinearInterpolateImageFunction< MovingImageType, double >       InterpolatorType;
	typedef itk::ImageRegistrationMethod< FixedImageType, MovingImageType >       RegistrationType;
	MetricType::Pointer         metric = MetricType::New();
	OptimizerType::Pointer      optimizer = OptimizerType::New();
	InterpolatorType::Pointer   interpolator = InterpolatorType::New();
	RegistrationType::Pointer   registration = RegistrationType::New();
	registration->SetMetric(metric);
	registration->SetOptimizer(optimizer);
	registration->SetInterpolator(interpolator);
	//  Software Guide : BeginLatex
	//
	//  The transform object is constructed below and passed to the registration
	//  method.
	//
	//  \index{itk::Versor\-Rigid3D\-Transform!New()}
	//  \index{itk::Versor\-Rigid3D\-Transform!Pointer}
	//  \index{itk::Registration\-Method!SetTransform()}
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	TransformType::Pointer  transform = TransformType::New();
	registration->SetTransform(transform);
	// Software Guide : EndCodeSnippet
	typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
	typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
	FixedImageReaderType::Pointer  fixedImageReader = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
	fixedImageReader->SetFileName("output_fixed.nrrd");
	movingImageReader->SetFileName("output_moving.nrrd");
	registration->SetFixedImage(fixedImageReader->GetOutput());
	registration->SetMovingImage(movingImageReader->GetOutput());
	fixedImageReader->Update();
	registration->SetFixedImageRegion(
		fixedImageReader->GetOutput()->GetBufferedRegion());
	//  Software Guide : BeginLatex
	//
	//  The input images are taken from readers. It is not necessary here to
	//  explicitly call \code{Update()} on the readers since the
	//  \doxygen{CenteredTransformInitializer} will do it as part of its
	//  computations. The following code instantiates the type of the
	//  initializer. This class is templated over the fixed and moving image type
	//  as well as the transform type. An initializer is then constructed by
	//  calling the \code{New()} method and assigning the result to a smart
	//  pointer.
	//
	// \index{itk::Centered\-Transform\-Initializer!Instantiation}
	// \index{itk::Centered\-Transform\-Initializer!New()}
	// \index{itk::Centered\-Transform\-Initializer!SmartPointer}
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	typedef itk::CenteredTransformInitializer< TransformType,
		FixedImageType,
		MovingImageType
	>  TransformInitializerType;
	TransformInitializerType::Pointer initializer =
		TransformInitializerType::New();
	// Software Guide : EndCodeSnippet
	//  Software Guide : BeginLatex
	//
	//  The initializer is now connected to the transform and to the fixed and
	//  moving images.
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	initializer->SetTransform(transform);
	initializer->SetFixedImage(fixedImageReader->GetOutput());
	initializer->SetMovingImage(movingImageReader->GetOutput());
	// Software Guide : EndCodeSnippet
	//  Software Guide : BeginLatex
	//
	//  The use of the geometrical centers is selected by calling
	//  \code{GeometryOn()} while the use of center of mass is selected by
	//  calling \code{MomentsOn()}.  Below we select the center of mass mode.
	//
	//  \index{Centered\-Transform\-Initializer!MomentsOn()}
	//  \index{Centered\-Transform\-Initializer!GeometryOn()}
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	initializer->MomentsOn();
	// Software Guide : EndCodeSnippet
	//  Software Guide : BeginLatex
	//
	//  Finally, the computation of the center and translation is triggered by
	//  the \code{InitializeTransform()} method. The resulting values will be
	//  passed directly to the transform.
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	initializer->InitializeTransform();
	// Software Guide : EndCodeSnippet
	//  Software Guide : BeginLatex
	//
	//  The rotation part of the transform is initialized using a
	//  \doxygen{Versor} which is simply a unit quaternion.  The
	//  \code{VersorType} can be obtained from the transform traits. The versor
	//  itself defines the type of the vector used to indicate the rotation axis.
	//  This trait can be extracted as \code{VectorType}. The following lines
	//  create a versor object and initialize its parameters by passing a
	//  rotation axis and an angle.
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	typedef TransformType::VersorType  VersorType;
	typedef VersorType::VectorType     VectorType;
	VersorType     rotation;
	VectorType     axis;
	axis[0] = 0.0;
	axis[1] = 0.0;
	axis[2] = 1.0;
	const double angle = 0;
	rotation.Set(axis, angle);
	transform->SetRotation(rotation);
	// Software Guide : EndCodeSnippet
	//  Software Guide : BeginLatex
	//
	//  We now pass the parameters of the current transform as the initial
	//  parameters to be used when the registration process starts.
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	registration->SetInitialTransformParameters(transform->GetParameters());
	// Software Guide : EndCodeSnippet
	typedef OptimizerType::ScalesType       OptimizerScalesType;
	OptimizerScalesType optimizerScales(transform->GetNumberOfParameters());
	const double translationScale = 1.0 / 1000.0;
	optimizerScales[0] = 1.0;
	optimizerScales[1] = 1.0;
	optimizerScales[2] = 1.0;
	optimizerScales[3] = translationScale;
	optimizerScales[4] = translationScale;
	optimizerScales[5] = translationScale;
	optimizer->SetScales(optimizerScales);
	optimizer->SetMaximumStepLength(0.2000);
	optimizer->SetMinimumStepLength(0.0001);
	optimizer->SetNumberOfIterations(200);
	// Create the Command observer and register it with the optimizer.
	//
	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
	optimizer->AddObserver(itk::IterationEvent(), observer);
	try
	{
		registration->Update();
		std::cout << "Optimizer stop condition: "
			<< registration->GetOptimizer()->GetStopConditionDescription()
			<< std::endl;
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
	}
	OptimizerType::ParametersType finalParameters =
		registration->GetLastTransformParameters();
	const double versorX = finalParameters[0];
	const double versorY = finalParameters[1];
	const double versorZ = finalParameters[2];
	const double finalTranslationX = finalParameters[3];
	const double finalTranslationY = finalParameters[4];
	const double finalTranslationZ = finalParameters[5];
	const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
	const double bestValue = optimizer->GetValue();
	// Print out results
	//
	std::cout << std::endl << std::endl;
	std::cout << "Result = " << std::endl;
	std::cout << " versor X      = " << versorX << std::endl;
	std::cout << " versor Y      = " << versorY << std::endl;
	std::cout << " versor Z      = " << versorZ << std::endl;
	std::cout << " Translation X = " << finalTranslationX << std::endl;
	std::cout << " Translation Y = " << finalTranslationY << std::endl;
	std::cout << " Translation Z = " << finalTranslationZ << std::endl;
	std::cout << " Iterations    = " << numberOfIterations << std::endl;
	std::cout << " Metric value  = " << bestValue << std::endl;

	const auto duration = std::chrono::duration_cast<ms>(clock::now() - before);

	std::cout << "It took " << duration.count() / 1000.0 << "ms" << std::endl;
}

void qSlicerVertRegModuleWidget::automaticRegistration() {
	int initialRotationX = 0;
	int initialRotationY = 0;
	int initialRotationZ = 20;

	int delta = 3;

	rotateAndProjectVolume(initialRotationX, initialRotationY, initialRotationZ);

	registerVolume2dEuler();


}

double qSlicerVertRegModuleWidget::registerVolume2dSimilarity() {
	constexpr unsigned int Dimension = 2;
	using PixelType = float;
	using FixedImageType = itk::Image< PixelType, Dimension >;
	using MovingImageType = itk::Image< PixelType, Dimension >;
	//  Software Guide : BeginLatex
	//
	//  The Transform class is instantiated using the code below. The only
	//  template parameter of this class is the representation type of the
	//  space coordinates.
	//
	//  \index{itk::Simularity2DTransform!Instantiation}
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	using TransformType = itk::Similarity2DTransform< double >;
	// Software Guide : EndCodeSnippet
	using OptimizerType = itk::RegularStepGradientDescentOptimizerv4<double>;
	using MetricType = itk::MeanSquaresImageToImageMetricv4< FixedImageType,
		MovingImageType >;
	using RegistrationType = itk::ImageRegistrationMethodv4< FixedImageType,
		MovingImageType,
		TransformType >;
	MetricType::Pointer         metric = MetricType::New();
	OptimizerType::Pointer      optimizer = OptimizerType::New();
	RegistrationType::Pointer   registration = RegistrationType::New();
	registration->SetMetric(metric);
	registration->SetOptimizer(optimizer);
	//  Software Guide : BeginLatex
	//
	//  As before, the transform object is constructed and initialized before it
	//  is passed to the registration filter.
	//
	//  \index{itk::Simularity2DTransform!Pointer}
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	TransformType::Pointer  transform = TransformType::New();
	// Software Guide : EndCodeSnippet
	using FixedImageReaderType = itk::ImageFileReader< FixedImageType  >;
	using MovingImageReaderType = itk::ImageFileReader< MovingImageType >;
	FixedImageReaderType::Pointer  fixedImageReader = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
	fixedImageReader->SetFileName("output_fixed.nrrd");
	movingImageReader->SetFileName("output_moving.nrrd");
	registration->SetFixedImage(fixedImageReader->GetOutput());
	registration->SetMovingImage(movingImageReader->GetOutput());
	//  Software Guide : BeginLatex
	//
	//  In this example, we again use the helper class
	//  \doxygen{CenteredTransformInitializer} to compute a reasonable
	//  value for the initial center of rotation and scaling along with
	//  an initial translation.
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	using TransformInitializerType = itk::CenteredTransformInitializer<
		TransformType,
		FixedImageType,
		MovingImageType >;
	TransformInitializerType::Pointer initializer
		= TransformInitializerType::New();
	initializer->SetTransform(transform);
	initializer->SetFixedImage(fixedImageReader->GetOutput());
	initializer->SetMovingImage(movingImageReader->GetOutput());
	initializer->MomentsOn();
	initializer->InitializeTransform();
	// Software Guide : EndCodeSnippet
	//  Software Guide : BeginLatex
	//
	//  The remaining parameters of the transform are initialized below.
	//
	//  \index{itk::Simularity2DTransform!SetScale()}
	//  \index{itk::Simularity2DTransform!SetAngle()}
	//
	//  Software Guide : EndLatex
	double initialScale = 1.0;
	/*if (argc > 7)
	{
		initialScale = std::stod(argv[7]);
	}*/
	double initialAngle = 0.0;
	/*if (argc > 8)
	{
		initialAngle = std::stod(argv[8]);
	}*/
	// Software Guide : BeginCodeSnippet
	transform->SetScale(initialScale);
	transform->SetAngle(initialAngle);
	// Software Guide : EndCodeSnippet
	//  Software Guide : BeginLatex
	//
	//  Now the initialized transform object will be set to the registration method,
	//  and its initial parameters are used to initialize the registration process.
	//
	//  Also, by calling the \code{InPlaceOn()} method, this initialized
	//  transform will be the output transform
	//  object or ``grafted'' to the output of the registration process.
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	registration->SetInitialTransform(transform);
	registration->InPlaceOn();
	// Software Guide : EndCodeSnippet
	//  Software Guide : BeginLatex
	//
	//  Keeping in mind that the scale of units in scaling, rotation and
	//  translation are quite different, we take advantage of the scaling
	//  functionality provided by the optimizers. We know that the first element
	//  of the parameters array corresponds to the scale factor, the second
	//  corresponds to the angle, third and fourth are the remaining
	//  translation. We use henceforth small factors in the scales
	//  associated with translations.
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	using OptimizerScalesType = OptimizerType::ScalesType;
	OptimizerScalesType optimizerScales(transform->GetNumberOfParameters());
	const double translationScale = 1.0 / 100.0;
	optimizerScales[0] = 10.0;
	optimizerScales[1] = 1.0;
	optimizerScales[2] = translationScale;
	optimizerScales[3] = translationScale;
	optimizer->SetScales(optimizerScales);
	// Software Guide : EndCodeSnippet
	//  Software Guide : BeginLatex
	//
	//  We also set the ordinary parameters of the optimization method. In this
	//  case we are using a
	//  \doxygen{RegularStepGradientDescentOptimizerv4}. Below we define the
	//  optimization parameters, i.e. initial learning rate (step length), minimal
	//  step length and number of iterations. The last two act as stopping criteria
	//  for the optimization.
	//
	//  Software Guide : EndLatex
	double steplength = 1.0;
	/*if (argc > 6)
	{
		steplength = std::stod(argv[6]);
	}*/
	// Software Guide : BeginCodeSnippet
	optimizer->SetLearningRate(steplength);
	optimizer->SetMinimumStepLength(0.0001);
	optimizer->SetNumberOfIterations(500);
	// Software Guide : EndCodeSnippet
	// Create the Command observer and register it with the optimizer.
	//
	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
	optimizer->AddObserver(itk::IterationEvent(), observer);
	// One level registration process without shrinking and smoothing.
	//
	constexpr unsigned int numberOfLevels = 1;
	RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
	shrinkFactorsPerLevel.SetSize(1);
	shrinkFactorsPerLevel[0] = 1;
	RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
	smoothingSigmasPerLevel.SetSize(1);
	smoothingSigmasPerLevel[0] = 0;
	registration->SetNumberOfLevels(numberOfLevels);
	registration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel);
	registration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel);
	try
	{
		registration->Update();
		std::cout << "Optimizer stop condition: "
			<< registration->GetOptimizer()->GetStopConditionDescription()
			<< std::endl;
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		//return EXIT_FAILURE;
	}
	TransformType::ParametersType finalParameters =
		transform->GetParameters();
	const double finalScale = finalParameters[0];
	const double finalAngle = finalParameters[1];
	const double finalTranslationX = finalParameters[2];
	const double finalTranslationY = finalParameters[3];
	const double rotationCenterX = registration->GetOutput()->Get()->GetFixedParameters()[0];
	const double rotationCenterY = registration->GetOutput()->Get()->GetFixedParameters()[1];
	const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
	const double bestValue = optimizer->GetValue();
	// Print out results
	//
	const double finalAngleInDegrees = finalAngle * 180.0 / itk::Math::pi;
	std::cout << std::endl;
	std::cout << "Result = " << std::endl;
	std::cout << " Scale           = " << finalScale << std::endl;
	std::cout << " Angle (radians) = " << finalAngle << std::endl;
	std::cout << " Angle (degrees) =  " << finalAngleInDegrees << std::endl;
	std::cout << " Translation X   = " << finalTranslationX << std::endl;
	std::cout << " Translation Y   = " << finalTranslationY << std::endl;
	std::cout << " Fixed Center X  = " << rotationCenterX << std::endl;
	std::cout << " Fixed Center Y  = " << rotationCenterY << std::endl;
	std::cout << " Iterations      = " << numberOfIterations << std::endl;
	std::cout << " Metric value    = " << bestValue << std::endl;

	return bestValue;
}

double qSlicerVertRegModuleWidget::registerVolume2dEuler() {
	std::cout << " Registration: Euler";

	constexpr unsigned int Dimension = 2;
	using PixelType = float;
	using FixedImageType = itk::Image< PixelType, Dimension >;
	using MovingImageType = itk::Image< PixelType, Dimension >;
	//  Software Guide : BeginLatex
	//
	//  The transform type is instantiated using the code below. The only
	//  template parameter of this class is the representation type of the
	//  space coordinates.
	//
	//  \index{itk::Euler2DTransform!Instantiation}
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	using TransformType = itk::Euler2DTransform< double >;
	// Software Guide : EndCodeSnippet
	using OptimizerType = itk::RegularStepGradientDescentOptimizerv4<double>;
	using MetricType = itk::MeanSquaresImageToImageMetricv4<
		FixedImageType,
		MovingImageType >;
	using RegistrationType = itk::ImageRegistrationMethodv4<
		FixedImageType,
		MovingImageType >;
	MetricType::Pointer         metric = MetricType::New();
	OptimizerType::Pointer      optimizer = OptimizerType::New();
	RegistrationType::Pointer   registration = RegistrationType::New();
	registration->SetMetric(metric);
	registration->SetOptimizer(optimizer);
	//  Software Guide : BeginLatex
	//
	//  Like the previous section, a direct initialization method is used here.
	//  The transform object is constructed below. This transform will
	//  be initialized, and its initial parameters will be considered as
	//  the parameters to be used when the registration process begins.
	//
	//  \index{itk::Euler2DTransform!Pointer}
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	TransformType::Pointer  transform = TransformType::New();
	// Software Guide : EndCodeSnippet
	using FixedImageReaderType = itk::ImageFileReader< FixedImageType  >;
	using MovingImageReaderType = itk::ImageFileReader< MovingImageType >;
	FixedImageReaderType::Pointer  fixedImageReader = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
	fixedImageReader->SetFileName("output_fixed.nrrd");
	movingImageReader->SetFileName("output_moving.nrrd");
	registration->SetFixedImage(fixedImageReader->GetOutput());
	registration->SetMovingImage(movingImageReader->GetOutput());
	//  Software Guide : BeginLatex
	//
	//  The input images are taken from readers. It is not necessary to
	//  explicitly call \code{Update()} on the readers since the
	//  CenteredTransformInitializer class will do it as part of its
	//  initialization. The following code instantiates the initializer. This
	//  class is templated over the fixed and moving images type as well as the
	//  transform type. An initializer is then constructed by calling the
	//  \code{New()} method and assigning the result to a
	//  \doxygen{SmartPointer}.
	//
	// \index{itk::Euler2DTransform!Instantiation}
	// \index{itk::Euler2DTransform!New()}
	// \index{itk::Euler2DTransform!SmartPointer}
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	using TransformInitializerType = itk::CenteredTransformInitializer<
		TransformType,
		FixedImageType,
		MovingImageType >;
	TransformInitializerType::Pointer initializer =
		TransformInitializerType::New();
	// Software Guide : EndCodeSnippet
	//  Software Guide : BeginLatex
	//
	//  The initializer is now connected to the transform and to the fixed and
	//  moving images.
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	initializer->SetTransform(transform);
	initializer->SetFixedImage(fixedImageReader->GetOutput());
	initializer->SetMovingImage(movingImageReader->GetOutput());
	// Software Guide : EndCodeSnippet
	//  Software Guide : BeginLatex
	//
	//  The use of the geometrical centers is selected by calling
	//  \code{GeometryOn()} while the use of center of mass is selected by
	//  calling \code{MomentsOn()}.  Below we select the center of mass mode.
	//
	//  \index{CenteredTransformInitializer!MomentsOn()}
	//  \index{CenteredTransformInitializer!GeometryOn()}
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	initializer->MomentsOn();
	// Software Guide : EndCodeSnippet
	//  Software Guide : BeginLatex
	//
	//  Finally, the computation of the center and translation is triggered by
	//  the \code{InitializeTransform()} method. The resulting values will be
	//  passed directly to the transform.
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	initializer->InitializeTransform();
	// Software Guide : EndCodeSnippet
	//  Software Guide : BeginLatex
	//
	//  The remaining parameters of the transform are initialized as before.
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	transform->SetAngle(0.0);
	// Software Guide : EndCodeSnippet
	//  Software Guide : BeginLatex
	//
	//  Now the initialized transform object will be set to the registration method,
	//  and the starting point of the registration is defined by its initial parameters.
	//
	//  If the \code{InPlaceOn()} method is called, this initialized transform will be the output transform
	//  object or ``grafted'' to the output. Otherwise, this ``InitialTransform'' will be deep-copied or
	//  ``cloned'' to the output.
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	registration->SetInitialTransform(transform);
	registration->InPlaceOn();
	// Software Guide : EndCodeSnippet
	using OptimizerScalesType = OptimizerType::ScalesType;
	OptimizerScalesType optimizerScales(transform->GetNumberOfParameters());
	const double translationScale = 1.0 / 1000.0;
	optimizerScales[0] = 1.0;
	optimizerScales[1] = translationScale;
	optimizerScales[2] = translationScale;
	optimizer->SetScales(optimizerScales);
	optimizer->SetLearningRate(0.1);
	optimizer->SetMinimumStepLength(0.001);
	optimizer->SetNumberOfIterations(200);
	// Create the Command observer and register it with the optimizer.
	//
	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
	optimizer->AddObserver(itk::IterationEvent(), observer);
	// One level registration process without shrinking and smoothing.
	//
	constexpr unsigned int numberOfLevels = 1;
	RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
	shrinkFactorsPerLevel.SetSize(1);
	shrinkFactorsPerLevel[0] = 1;
	RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
	smoothingSigmasPerLevel.SetSize(1);
	smoothingSigmasPerLevel[0] = 0;
	registration->SetNumberOfLevels(numberOfLevels);
	registration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel);
	registration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel);
	try
	{
		registration->Update();
		std::cout << "Optimizer stop condition: "
			<< registration->GetOptimizer()->GetStopConditionDescription()
			<< std::endl;
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		//return EXIT_FAILURE;
	}
	//  Software Guide : BeginLatex
	//
	//  Since the registration filter has \code{InPlace} set, the transform object
	//  is grafted to the output and is updated by the registration method.
	//
	//  Software Guide : EndLatex
	TransformType::ParametersType finalParameters = transform->GetParameters();
	const double finalAngle = finalParameters[0];
	const double finalTranslationX = finalParameters[1];
	const double finalTranslationY = finalParameters[2];
	const double rotationCenterX = registration->GetOutput()->Get()->GetFixedParameters()[0];
	const double rotationCenterY = registration->GetOutput()->Get()->GetFixedParameters()[1];
	const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
	const double bestValue = optimizer->GetValue();
	// Print out results
	//
	const double finalAngleInDegrees = finalAngle * 180.0 / itk::Math::pi;
	std::cout << "Result = " << std::endl;
	std::cout << " Angle (radians) " << finalAngle << std::endl;
	std::cout << " Angle (degrees)  " << finalAngleInDegrees << std::endl;
	std::cout << " Translation X  = " << finalTranslationX << std::endl;
	std::cout << " Translation Y  = " << finalTranslationY << std::endl;
	std::cout << " Fixed Center X = " << rotationCenterX << std::endl;
	std::cout << " Fixed Center Y = " << rotationCenterY << std::endl;
	std::cout << " Iterations     = " << numberOfIterations << std::endl;
	std::cout << " Metric value   = " << bestValue << std::endl;
	//  Software Guide : BeginLatex
	//
	//  Let's execute this example over some of the images provided in
	//  \code{Examples/Data}, for example:
	//
	//  \begin{itemize}
	//  \item \code{BrainProtonDensitySliceBorder20.png}
	//  \item \code{BrainProtonDensitySliceR10X13Y17.png}
	//  \end{itemize}
	//
	//  The second image is the result of intentionally rotating the first
	//  image by $10$ degrees around the geometric center and shifting
	//  it $13mm$ in $X$ and $17mm$ in $Y$. Both images have
	//  unit-spacing and are shown in Figure
	//  \ref{fig:FixedMovingImageRegistration5}. The registration takes
	//  $21$ iterations and produces:
	//
	//  \begin{center}
	//  \begin{verbatim}
	//  [ 0.174527, 12.4528, 16.0766]
	//  \end{verbatim}
	//  \end{center}
	//
	//  These parameters are interpreted as
	//
	//  \begin{itemize}
	//  \item Angle         =                  $0.174527$     radians
	//  \item Translation   = $( 12.4528, 16.0766 )$ millimeters
	//  \end{itemize}
	//
	//  Note that the reported translation is not the translation of $(13,17)$
	//  that might be expected. The reason is that we used the center of
	//  mass $( 111.204, 131.591 )$  for the fixed center, while the input was rotated
	//  about the geometric center $( 110.5, 128.5 )$.  It is more illustrative in
	//  this case to take a look at the actual rotation matrix and offset
	//  resulting from the five parameters.
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	TransformType::MatrixType matrix = transform->GetMatrix();
	TransformType::OffsetType offset = transform->GetOffset();
	std::cout << "Matrix = " << std::endl << matrix << std::endl;
	std::cout << "Offset = " << std::endl << offset << std::endl;

	using ResampleFilterType = itk::ResampleImageFilter<
		MovingImageType,
		FixedImageType >;
	ResampleFilterType::Pointer resample = ResampleFilterType::New();
	resample->SetTransform(transform);
	resample->SetInput(movingImageReader->GetOutput());
	FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();
	resample->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
	resample->SetOutputOrigin(fixedImage->GetOrigin());
	resample->SetOutputSpacing(fixedImage->GetSpacing());
	resample->SetOutputDirection(fixedImage->GetDirection());
	resample->SetDefaultPixelValue(100);
	using OutputPixelType = unsigned char;
	using OutputImageType = itk::Image< OutputPixelType, Dimension >;
	using CastFilterType = itk::CastImageFilter<
		FixedImageType,
		OutputImageType >;
	using WriterType = itk::ImageFileWriter< OutputImageType >;
	WriterType::Pointer      writer = WriterType::New();
	CastFilterType::Pointer  caster = CastFilterType::New();
	writer->SetFileName("output_euler.nrrd");
	caster->SetInput(resample->GetOutput());
	writer->SetInput(caster->GetOutput());
	writer->Update();

	// Now compute the difference between the images
	// before and after registration.
	//
	using DifferenceImageType = itk::Image< float, Dimension >;
	using DifferenceFilterType = itk::SubtractImageFilter<
		FixedImageType,
		FixedImageType,
		DifferenceImageType >;
	DifferenceFilterType::Pointer difference = DifferenceFilterType::New();
	using OutputPixelType = unsigned char;
	using OutputImageType = itk::Image< OutputPixelType, Dimension >;
	using RescalerType = itk::RescaleIntensityImageFilter<
		DifferenceImageType,
		OutputImageType >;
	RescalerType::Pointer intensityRescaler = RescalerType::New();
	intensityRescaler->SetOutputMinimum(0);
	intensityRescaler->SetOutputMaximum(255);
	difference->SetInput1(fixedImageReader->GetOutput());
	difference->SetInput2(resample->GetOutput());
	resample->SetDefaultPixelValue(1);
	intensityRescaler->SetInput(difference->GetOutput());
	using WriterType = itk::ImageFileWriter< OutputImageType >;
	WriterType::Pointer      writer2 = WriterType::New();
	writer2->SetInput(intensityRescaler->GetOutput());
	try
	{
		// Compute the difference image between the
		// fixed and moving image after registration.
		//if (argc > 5)
		//{
			writer2->SetFileName("output_euler_before.nrrd");
			writer2->Update();
		//}
		// Compute the difference image between the
		// fixed and resampled moving image after registration.
		TransformType::Pointer identityTransform = TransformType::New();
		identityTransform->SetIdentity();
		resample->SetTransform(identityTransform);
		//if (argc > 4)
		//{
			writer2->SetFileName("output_euler_after.nrrd");
			writer2->Update();
		//}
	}
	catch (itk::ExceptionObject & excp)
	{
		std::cerr << "Error while writing difference images" << std::endl;
		std::cerr << excp << std::endl;
		//return EXIT_FAILURE;
	}
	//return EXIT_SUCCESS;

	return bestValue;
}

void qSlicerVertRegModuleWidget::registerVolumeOld()
{
	Q_D(qSlicerVertRegModuleWidget);

	using clock = std::chrono::system_clock;
	using ms = std::chrono::milliseconds;

	const auto before = clock::now();

	qDebug() << "function::registerVolume()";

	typedef itk::Image<double, 2>                 FixedImageType;
	typedef itk::Image<double, 2>                 MovingImageType;
	typedef itk::ImageFileReader<FixedImageType>  FixedImageReaderType;
	typedef itk::ImageFileReader<MovingImageType> MovingImageReaderType;
	typedef itk::Euler2DTransform< double >       TransformType;
	typedef itk::ExhaustiveOptimizerv4< double >  OptimizerType;
	typedef itk::MeanSquaresImageToImageMetricv4< FixedImageType,
		MovingImageType >
		MetricType;
	typedef itk::CenteredTransformInitializer< TransformType,
		FixedImageType, MovingImageType >
		TransformInitializerType;
	typedef itk::ImageRegistrationMethodv4< FixedImageType,
		MovingImageType, TransformType >
		RegistrationType;

	FixedImageType::Pointer           fixedImage = FixedImageType::New();
	MovingImageType::Pointer          movingImage = MovingImageType::New();
	TransformType::Pointer            transform = TransformType::New();
	MetricType::Pointer               metric = MetricType::New();
	OptimizerType::Pointer            optimizer = OptimizerType::New();
	RegistrationType::Pointer         registration = RegistrationType::New();
	TransformInitializerType::Pointer initializer =
		TransformInitializerType::New();

	vtkSmartPointer<vtkImageData> inputVTK = vtkSmartPointer<vtkImageData>::New();
	inputVTK = d->InputVolumeNode->GetImageData();
	vtkSmartPointer<vtkImageCast> castFilter =
	vtkSmartPointer<vtkImageCast>::New();
	castFilter->SetOutputScalarTypeToDouble();
	castFilter->SetInputData(inputVTK);
	castFilter->Update();
	typedef itk::VTKImageToImageFilter< FixedImageType > VTKtoITKFilterType;
	VTKtoITKFilterType::Pointer vtkToItkFilter = VTKtoITKFilterType::New();
	vtkToItkFilter->SetInput(castFilter->GetOutput());
	vtkToItkFilter->Update();
	fixedImage = vtkToItkFilter->GetOutput();

	vtkSmartPointer<vtkImageData> inputVTK2 = vtkSmartPointer<vtkImageData>::New();
	inputVTK2 = d->OutputVolumeNode->GetImageData();
	vtkSmartPointer<vtkImageCast> castFilter2 =
		vtkSmartPointer<vtkImageCast>::New();
	castFilter2->SetOutputScalarTypeToDouble();
	castFilter2->SetInputData(inputVTK2);
	castFilter2->Update();
	typedef itk::VTKImageToImageFilter< FixedImageType > VTKtoITKFilterType2;
	VTKtoITKFilterType2::Pointer vtkToItkFilter2 = VTKtoITKFilterType2::New();
	vtkToItkFilter2->SetInput(castFilter2->GetOutput());
	vtkToItkFilter2->Update();
	movingImage = vtkToItkFilter->GetOutput();

	unsigned int angles = 12;
	OptimizerType::StepsType steps(transform->GetNumberOfParameters());
	steps[0] = int(angles / 2);
	steps[1] = 0;
	steps[2] = 0;
	optimizer->SetNumberOfSteps(steps);

	OptimizerType::ScalesType scales(transform->GetNumberOfParameters());
	scales[0] = 2.0 * vnl_math::pi / angles;
	scales[1] = 1.0;
	scales[2] = 1.0;

	optimizer->SetScales(scales);

	initializer->SetTransform(transform);
	initializer->SetFixedImage(fixedImage);
	initializer->SetMovingImage(movingImage);
	initializer->InitializeTransform();

	// Initialize registration
	registration->SetMetric(metric);
	registration->SetOptimizer(optimizer);
	registration->SetFixedImage(fixedImage);
	registration->SetMovingImage(movingImage);
	registration->SetInitialTransform(transform);

	try
	{
		registration->Update();
		optimizer->Print(std::cout);
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
	}

	const auto duration = std::chrono::duration_cast<ms>(clock::now() - before);

	std::cout << "It took " << duration.count() / 1000.0 << "ms" << std::endl;
}

void qSlicerVertRegModuleWidget::rotateAndProjectVolume0() {
	rotateAndProjectVolume(0, 0, 0);
}

void qSlicerVertRegModuleWidget::rotateAndProjectVolume5() {
	rotateAndProjectVolume(0, 0, 5);
}

void qSlicerVertRegModuleWidget::rotateAndProjectVolume15() {
	rotateAndProjectVolume(0, 0, 15);
}

void qSlicerVertRegModuleWidget::rotateAndProjectVolume25() {
	rotateAndProjectVolume(0, 0, 25);
}

void qSlicerVertRegModuleWidget::rotateAndProjectVolume(int x, int y, int z)
{
	Q_D(qSlicerVertRegModuleWidget);

	using clock = std::chrono::system_clock;
	using ms = std::chrono::milliseconds;

	const auto before = clock::now();

	qDebug() << "function::rotateAndProjectVolume()";

	//int x = 0;
	//int y = 0;
	//int z = 20;

	const int dim = 3;
	typedef short PType;

	typedef itk::Image< PType, dim > InputImageType;
	typedef itk::Image< PType, dim - 1 > OutputImageType;

	vtkSmartPointer<vtkImageData> inputVTK = vtkSmartPointer<vtkImageData>::New();
	inputVTK = d->InputVolumeNode->GetImageData();

	//vtkNew<vtkMatrix4x4> matr;
	//d->InputVolumeNode->GetIJKToRASMatrix(matr);
	//d->OutputVolumeNode->SetIJKToRASMatrix(matr);

	double bounds29[6];
	inputVTK->GetBounds(bounds29);

	// Rotate about the center of the image
	vtkSmartPointer<vtkTransform> transform29 =
		vtkSmartPointer<vtkTransform>::New();

	// Compute the center of the image
	double center29[3];
	center29[0] = (bounds29[1] + bounds29[0]) / 2.0;
	center29[1] = (bounds29[3] + bounds29[2]) / 2.0;
	center29[2] = (bounds29[5] + bounds29[4]) / 2.0;

	// Rotate about the center
	transform29->Translate(center29[0], center29[1], center29[2]);
	transform29->RotateWXYZ((180+z), 0, 0, 1);
	transform29->Translate(-center29[0], -center29[1], -center29[2]);

	// Reslice does all of the work
	vtkSmartPointer<vtkImageReslice> reslice29 =
		vtkSmartPointer<vtkImageReslice>::New();
	reslice29->SetInputData(inputVTK);
	reslice29->SetResliceTransform(transform29);
	reslice29->SetInterpolationModeToCubic();
	reslice29->SetOutputSpacing(
		inputVTK->GetSpacing()[0],
		inputVTK->GetSpacing()[1],
		inputVTK->GetSpacing()[2]);
	reslice29->SetOutputOrigin(
		inputVTK->GetOrigin()[0],
		inputVTK->GetOrigin()[1],
		inputVTK->GetOrigin()[2]);
	reslice29->SetOutputExtent(inputVTK->GetExtent()); // Use a larger extent than the original image's to prevent clipping
	reslice29->Update();

	qDebug() << "function::rotateAndProjectVolume1()";

	typedef itk::VTKImageToImageFilter< InputImageType > VTKtoITKFilterType;
	VTKtoITKFilterType::Pointer vtkToItkFilter = VTKtoITKFilterType::New();
	//vtkToItkFilter->SetInput(inputVTK);
	vtkToItkFilter->SetInput(reslice29->GetOutput());
	vtkToItkFilter->Update();
	InputImageType::Pointer image = vtkToItkFilter->GetOutput();

	qDebug() << "function::rotateAndProjectVolume2()";

	/*typedef itk::CenteredAffineTransform< double, dim > TransformType;
	typedef itk::ResampleImageFilter<InputImageType, InputImageType> RotateFilterType;
	RotateFilterType::Pointer rotatefilter = RotateFilterType::New();
	TransformType::Pointer transform = TransformType::New();

	TransformType::InputPointType origin;
	InputImageType::SizeType size;


	rotatefilter->SetTransform(transform);
	size = image->GetLargestPossibleRegion().GetSize();

	rotatefilter->SetInput(image);
	TransformType::InputPointType center;

	origin = image->GetOrigin();
	center[0] = origin[0] + size[0] / 2.0;
	center[1] = origin[1] + size[1] / 2.0;
	center[2] = origin[2] + size[2] / 2.0;

	transform->SetCenter(center);
	const double degreesToRadians = atan(1.0) / 45.0;
	double angle = x * degreesToRadians;
	TransformType::OutputVectorType axis;
	axis[0] = 1;
	axis[1] = 0;
	axis[2] = 0;

	transform->Rotate3D(axis, angle, false);

	angle = y * degreesToRadians;
	axis[0] = 0;
	axis[1] = 1;
	axis[2] = 0;

	transform->Rotate3D(axis, angle, false);

	angle = z * degreesToRadians;
	axis[0] = 0;
	axis[1] = 0;
	axis[2] = 1 + (180 * degreesToRadians);

	transform->Rotate3D(axis, angle, false);

	rotatefilter->SetDefaultPixelValue(100);
	rotatefilter->SetSize(size);
	rotatefilter->Update();*/

	using FilterType = itk::ThresholdImageFilter< InputImageType >;
	FilterType::Pointer thresholdfilter = FilterType::New();
	//thresholdfilter->SetInput(rotatefilter->GetOutput());
	thresholdfilter->SetInput(vtkToItkFilter->GetOutput());
	thresholdfilter->ThresholdOutside(95, 500);
	thresholdfilter->SetOutsideValue(0);
	thresholdfilter->Update();

	typedef itk::MeanProjectionImageFilter< InputImageType, OutputImageType > ProjectionFilterType;
	ProjectionFilterType::Pointer projfilter = ProjectionFilterType::New();
	projfilter->SetInput(thresholdfilter->GetOutput());
	projfilter->SetProjectionDimension(1);
	projfilter->Update();

	//itk::SimpleFilterWatcher watcher(projfilter, "filter");
	qDebug() << "function::rotateAndProjectVolume3()";


	//OutputImageType::Pointer outputITK = projfilter->GetOutput();

	qDebug() << "function::projectVolumeFilter4()";

	using RescaleFilterType = itk::RescaleIntensityImageFilter<
		OutputImageType, OutputImageType >;
	RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255); //TODO USED TO BE 255
	rescaler->SetInput(projfilter->GetOutput());
	rescaler->Update();

	OutputImageType::Pointer outputITK = rescaler->GetOutput();

	if (debug) {
		typedef itk::ImageFileWriter< OutputImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetInput(outputITK);
		writer->SetFileName("output_moving.nrrd");
		writer->Update();
	}

	qDebug() << "function::projectVolumeFilter5()";

	vtkSmartPointer<vtkImageData> outputVTK = vtkSmartPointer<vtkImageData>::New();

	qDebug() << "function::projectVolumeFilter6()";

	if (debug) {

		vtkSmartPointer<vtkNrrdReader> writer = vtkSmartPointer<vtkNrrdReader>::New();
		writer->SetInputData(outputVTK);
		writer->SetFileName("output_moving.nrrd");
		writer->Update();
		outputVTK = writer->GetOutput();
	}

	d->OutputVolumeNode->SetAndObserveImageData(outputVTK);

	const auto duration = std::chrono::duration_cast<ms>(clock::now() - before);

	std::cout << "It took " << duration.count() / 1000.0 << "ms" << std::endl;
}
//-----------------------------------------------------------------------------

void qSlicerVertRegModuleWidget::projectVolumeDRR0() {
	projectVolumeDRR3(0, 0, 0);
}

void qSlicerVertRegModuleWidget::projectVolumeDRR5() {
	projectVolumeDRR3(0, 0, 5);
}

void qSlicerVertRegModuleWidget::projectVolumeDRR15() {
	projectVolumeDRR3(0, 0, 15);
}

void qSlicerVertRegModuleWidget::projectVolumeDRR25() {
	projectVolumeDRR3(0, 0, 25);
}

void qSlicerVertRegModuleWidget::projectVolumeDRR3(int x, int y, int z)
{
	Q_D(qSlicerVertRegModuleWidget);

	using clock = std::chrono::system_clock;
	using ms = std::chrono::milliseconds;

	const auto before = clock::now();

	qDebug() << "function::projectVolumeDRR()";

	const int dim = 3;

	typedef short PType;
	typedef itk::Image< PType, dim > InputImageType;
	typedef itk::Image< PType, 2 > OutputImageType;


	vtkSmartPointer<vtkImageData> inputVTK = vtkSmartPointer<vtkImageData>::New();
	inputVTK = d->InputVolumeNode->GetImageData();

	//vtkNew<vtkMatrix4x4> matr;
	//d->InputVolumeNode->GetIJKToRASMatrix(matr);
	//d->OutputVolumeNode->SetIJKToRASMatrix(matr);

	double bounds29[6];
	inputVTK->GetBounds(bounds29);

	// Rotate about the center of the image
	vtkSmartPointer<vtkTransform> transform29 =
		vtkSmartPointer<vtkTransform>::New();

	// Compute the center of the image
	double center29[3];
	center29[0] = (bounds29[1] + bounds29[0]) / 2.0;
	center29[1] = (bounds29[3] + bounds29[2]) / 2.0;
	center29[2] = (bounds29[5] + bounds29[4]) / 2.0;

	// Rotate about the center
	transform29->Translate(center29[0], center29[1], center29[2]);
	transform29->RotateWXYZ((180 + z), 0, 0, 1);
	transform29->Translate(-center29[0], -center29[1], -center29[2]);

	// Reslice does all of the work
	vtkSmartPointer<vtkImageReslice> reslice29 =
		vtkSmartPointer<vtkImageReslice>::New();
	reslice29->SetInputData(inputVTK);
	reslice29->SetResliceTransform(transform29);
	reslice29->SetInterpolationModeToCubic();
	reslice29->SetOutputSpacing(
		inputVTK->GetSpacing()[0],
		inputVTK->GetSpacing()[1],
		inputVTK->GetSpacing()[2]);
	reslice29->SetOutputOrigin(
		inputVTK->GetOrigin()[0],
		inputVTK->GetOrigin()[1],
		inputVTK->GetOrigin()[2]);
	reslice29->SetOutputExtent(inputVTK->GetExtent()); // Use a larger extent than the original image's to prevent clipping
	reslice29->Update();

	// Convert VTK -> ITK
	typedef itk::VTKImageToImageFilter< InputImageType > VTKtoITKFilterType;
	VTKtoITKFilterType::Pointer vtkToItkFilter = VTKtoITKFilterType::New();
	//vtkToItkFilter->SetInput(inputVTK);
	vtkToItkFilter->SetInput(reslice29->GetOutput());
	vtkToItkFilter->Update();
	InputImageType::Pointer image = vtkToItkFilter->GetOutput();

	// Software Guide : BeginLatex
	//
	// Creation of a \code{ResampleImageFilter} enables coordinates for
	// each of the pixels in the DRR image to be generated. These
	// coordinates are used by the \code{RayCastInterpolateImageFunction}
	// to determine the equation of each corresponding ray which is cast
	// through the input volume.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	using FilterType = itk::ResampleImageFilter<InputImageType, InputImageType >;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput(image);
	filter->SetDefaultPixelValue(0);
	//qDebug() << "function::projectVolume2()::ResampleImageFilter";

	// Software Guide : EndCodeSnippet
	// Software Guide : BeginLatex
	//
	// An Euler transformation is defined to position the input volume.
	// The \code{ResampleImageFilter} uses this transform to position the
	// output DRR image for the desired view.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	using TransformType = itk::CenteredEuler3DTransform< double >;
	TransformType::Pointer transform = TransformType::New();
	transform->SetComputeZYX(true);
	TransformType::OutputVectorType translation;
	translation[0] = 0;// tx; left -right Increasing this moves the model right
	translation[1] = 200;// ty; THIS WAS ORIGINALLY 0 Increasing this zooms out further
	translation[2] = 0;// tz; up -down Increasing this moves the model downwards
					   // constant for converting degrees into radians
	const double dtr = (std::atan(1.0) * 4.0) / 180.0;
	transform->SetTranslation(translation);
	//transform->SetRotation(dtr*rx, dtr*ry, dtr*rz);
	//transform->SetRotation(dtr*0., dtr*0., dtr*0.);
	transform->SetRotation(dtr*90., dtr*0., dtr*0.); //TODO change 3rd coordinate to 20
													 //transform->SetRotation(dtr*(90.+x), dtr*y, dtr*z);
													 //forward-back, clockwise-ccw, 


	InputImageType::PointType   imOrigin = image->GetOrigin();
	InputImageType::SpacingType imRes = image->GetSpacing();
	using InputImageRegionType = InputImageType::RegionType;
	using InputImageSizeType = InputImageRegionType::SizeType;
	InputImageRegionType imRegion = image->GetBufferedRegion();
	InputImageSizeType   imSize = imRegion.GetSize();
	imOrigin[0] += imRes[0] * static_cast<double>(imSize[0]) / 2.0;
	imOrigin[1] += imRes[1] * static_cast<double>(imSize[1]) / 2.0;
	imOrigin[2] += imRes[2] * static_cast<double>(imSize[2]) / 2.0;
	TransformType::InputPointType center;
	/*center[0] = cx + imOrigin[0];
	center[1] = cy + imOrigin[1];
	center[2] = cz + imOrigin[2];
	transform->SetCenter(center);*/
	center[0] = 0 + imOrigin[0];
	center[1] = 0 + imOrigin[1];
	center[2] = 0 + imOrigin[2];
	transform->SetCenter(center);
	// Software Guide : EndCodeSnippet
	// Software Guide : BeginLatex
	//
	// The \code{RayCastInterpolateImageFunction} is instantiated and passed the transform
	// object. The \code{RayCastInterpolateImageFunction} uses this
	// transform to reposition the x-ray source such that the DRR image
	// and x-ray source move as one around the input volume. This coupling
	// mimics the rigid geometry of the x-ray gantry.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	using InterpolatorType =
		itk::RayCastInterpolateImageFunction<InputImageType, double>;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	interpolator->SetTransform(transform);
	// Software Guide : EndCodeSnippet
	// Software Guide : BeginLatex
	//
	// We can then specify a threshold above which the volume's
	// intensities will be integrated.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet

	interpolator->SetThreshold(0);//threshold);

								  // Software Guide : EndCodeSnippet
								  // Software Guide : BeginLatex
								  //
								  // The ray-cast interpolator needs to know the initial position of the
								  // ray source or focal point. In this example we place the input
								  // volume at the origin and halfway between the ray source and the
								  // screen. The distance between the ray source and the screen
								  // is the "source to image distance" \code{sid} and is specified by
								  // the user.
								  //
								  // Software Guide : EndLatex
								  // Software Guide : BeginCodeSnippet
	InterpolatorType::InputPointType focalpoint;
	focalpoint[0] = imOrigin[0] ;
	focalpoint[1] = imOrigin[1] ;
	focalpoint[2] = imOrigin[2] - 400 / 2 ;//- sid / 2.;
	interpolator->SetFocalPoint(focalpoint);
	// Software Guide : EndCodeSnippet
	// Software Guide : BeginLatex
	//
	// Having initialised the interpolator we pass the object to the
	// resample filter.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	//interpolator->Print(std::cout);
	filter->SetInterpolator(interpolator);
	filter->SetTransform(transform);
	// Software Guide : EndCodeSnippet
	// Software Guide : BeginLatex
	//
	// The size and resolution of the output DRR image is specified via the
	// resample filter.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	// setup the scene
	InputImageType::SizeType   size;
	size[0] = 512;//dx;  // number of pixels along X of the 2D DRR image
	size[1] = 494;// dy;  // number of pixels along Y of the 2D DRR image
	size[2] = 1;   // only one slice
	filter->SetSize(size);
	InputImageType::SpacingType spacing;
	spacing[0] = 1;// sx;  // pixel spacing along X of the 2D DRR image [mm]
	spacing[1] = 1;// sy;  // pixel spacing along Y of the 2D DRR image [mm]
	spacing[2] = 1; // slice thickness of the 2D DRR image [mm]
	filter->SetOutputSpacing(spacing);
	// Software Guide : EndCodeSnippet

	// Software Guide : BeginLatex
	//
	// In addition the position of the DRR is specified. The default
	// position of the input volume, prior to its transformation is
	// half-way between the ray source and screen and unless specified
	// otherwise the normal from the "screen" to the ray source passes
	// directly through the centre of the DRR.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	double origin[dim];
	origin[0] = imOrigin[0] + 0 - 1.*((double)512 - 1.) / 2.;//o2Dx - sx*((double)dx - 1.) / 2.;
	origin[1] = imOrigin[1] + 0 - 1.*((double)494 - 1.) / 2.;//o2Dy - sy*((double)dy - 1.) / 2.;
	origin[2] = imOrigin[2] + 400 / 2.;//sid / 2.;
	filter->SetOutputOrigin(origin);
	// Software Guide : EndCodeSnippet
	// create writer
	//qDebug() << "function::projectVolume2()::RayCastImageFilter";

	// Software Guide : BeginLatex
	//
	// The output of the resample filter can then be passed to a writer to
	// save the DRR image to a file.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	using RescaleFilterType = itk::RescaleIntensityImageFilter<
		InputImageType, InputImageType >;
	RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255); //TODO USED TO BE 255
	rescaler->SetInput(filter->GetOutput());
	rescaler->Update();

	typedef itk::ImageFileWriter< InputImageType >  WriterType;
	WriterType::Pointer writer9 = WriterType::New();
	writer9->SetFileName("output_fixed.nrrd");
	writer9->SetInput(rescaler->GetOutput());
	writer9->Update();

	typedef itk::ExtractImageFilter< InputImageType, OutputImageType > SliceType3;
	SliceType3::Pointer slice3 = SliceType3::New();
	slice3->InPlaceOn();
	slice3->SetDirectionCollapseToSubmatrix();
	InputImageType::RegionType inputRegion3 =
		rescaler->GetOutput()->GetLargestPossibleRegion();
	InputImageType::SizeType size3 = inputRegion3.GetSize();
	size3[2] = 0;
	InputImageType::IndexType start3 = inputRegion3.GetIndex();
	const unsigned int sliceNumber3 = 0;
	start3[2] = sliceNumber3;
	InputImageType::RegionType desiredRegion3;
	desiredRegion3.SetSize(size3);
	desiredRegion3.SetIndex(start3);
	slice3->SetExtractionRegion(desiredRegion3);
	//qDebug() << "function::projectVolume2()::" << size3[0] << " " << size3[1] << " ~ " << start3[0] << " " << start3[1];
	slice3->SetInput(rescaler->GetOutput());
	slice3->Update();

	/*WriterType::Pointer writer4 = WriterType::New();
	writer4->SetFileName("output4.jpeg");

	using ExtractFilterType = itk::ExtractImageFilter< OutputImageType, OutputImageType >;
	ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
	extractFilter->SetDirectionCollapseToSubmatrix();

	// set up the extraction region [one slice]
	const OutputImageType * inputImage = rescaler->GetOutput();
	OutputImageType::RegionType inputRegion = inputImage->GetBufferedRegion();
	OutputImageType::SizeType size1 = inputRegion.GetSize();
	size1[2] = 1; // we extract along z direction
	OutputImageType::IndexType start = inputRegion.GetIndex();
	const unsigned int sliceNumber = 0;//atoi(argv[3]);
	start[2] = sliceNumber;
	OutputImageType::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);

	extractFilter->SetExtractionRegion(desiredRegion);
	using PasteFilterType = itk::PasteImageFilter< OutputImageType >;
	PasteFilterType::Pointer pasteFilter = PasteFilterType::New();
	using MedianFilterType = itk::MedianImageFilter< OutputImageType,
	OutputImageType >;
	MedianFilterType::Pointer medianFilter = MedianFilterType::New();
	extractFilter->SetInput(inputImage);
	medianFilter->SetInput(extractFilter->GetOutput());
	pasteFilter->SetSourceImage(medianFilter->GetOutput());
	pasteFilter->SetDestinationImage(inputImage);
	pasteFilter->SetDestinationIndex(start);

	OutputImageType::SizeType indexRadius;
	indexRadius[0] = 1; // radius along x
	indexRadius[1] = 1; // radius along y
	indexRadius[2] = 0; // radius along z
	medianFilter->SetRadius(indexRadius);
	medianFilter->UpdateLargestPossibleRegion();
	const OutputImageType * medianImage = medianFilter->GetOutput();
	pasteFilter->SetSourceRegion(medianImage->GetBufferedRegion());
	writer4->SetInput(pasteFilter->GetOutput());
	writer4->Update();*/

	/*OutputImageType::IndexType start1;
	start1[0] = 0;  // first index on X
	start1[1] = 0;  // first index on Y
	start1[2] = 0;  // first index on Z
	OutputImageType::SizeType  size1;
	size1[0] = 490;  // size along X
	size1[1] = 490;  // size along Y
	size1[2] = 1;  // size along Z
	OutputImageType::RegionType region;
	region.SetSize(size1);
	region.SetIndex(start1);

	typedef itk::ExtractImageFilter< OutputImageType, OutputImageType2 > ExtractImageFilterType;
	ExtractImageFilterType::Pointer extract = ExtractImageFilterType::New();
	extract->SetInput(rescaler->GetOutput());
	extract->SetDirectionCollapseToIdentity();
	extract->SetExtractionRegion(region);
	extract->Update();
	using WriterType4 = itk::ImageFileWriter< OutputImageType2 >;
	WriterType4::Pointer writer4 = WriterType4::New();
	writer4->SetFileName("output4.jpeg");
	writer4->SetInput(extract->GetOutput());*/

	using WriterType11 = itk::ImageFileWriter< OutputImageType >;
	WriterType11::Pointer writer11 = WriterType11::New();
	writer11->SetFileName("output_fixed2.nrrd");
	writer11->SetInput(slice3->GetOutput());
	writer11->Update();


	using ITKtoVTKFilterType = itk::ImageToVTKImageFilter< OutputImageType >;
	ITKtoVTKFilterType::Pointer itktovtkfilter = ITKtoVTKFilterType::New();
	itktovtkfilter->SetInput(slice3->GetOutput());
	itktovtkfilter->Update();
	vtkImageData * outputVTK = itktovtkfilter->GetOutput();

	/*vtkSmartPointer<vtkImageReslice> reslice =
	vtkSmartPointer<vtkImageReslice>::New();
	reslice->SetOutputExtent(0, 9, 0, 100, 0, 0);
	reslice->SetInputData(outputVTK);
	reslice->Update();*/

	/*vtkSmartPointer<vtkImageWriter> writer5 =
	vtkSmartPointer<vtkImageWriter>::New();
	writer5->SetInputData(outputVTK);
	writer5->SetFileName("output3.vtk");
	writer5->Write();*/

	vtkSmartPointer<vtkImageCast> castFilter =
		vtkSmartPointer<vtkImageCast>::New();
	castFilter->SetOutputScalarTypeToUnsignedChar();
	castFilter->SetInputData(outputVTK);
	castFilter->Update();

	d->OutputVolumeNode->SetAndObserveImageData(castFilter->GetOutput());

	const auto duration = std::chrono::duration_cast<ms>(clock::now() - before);

	std::cout << "It took " << duration.count() / 1000.0 << "s" << std::endl;
}

void qSlicerVertRegModuleWidget::projectVolumeDRR3backup(int x, int y, int z)
{
	Q_D(qSlicerVertRegModuleWidget);

	using clock = std::chrono::system_clock;
	using ms = std::chrono::milliseconds;

	const auto before = clock::now();

	qDebug() << "function::projectVolumeDRR()";

	const int dim = 3;

	typedef short PType;
	typedef itk::Image< PType, dim > InputImageType;
	typedef itk::Image< PType, 2 > OutputImageType;


	vtkSmartPointer<vtkImageData> inputVTK = vtkSmartPointer<vtkImageData>::New();
	inputVTK = d->InputVolumeNode->GetImageData();

	//vtkNew<vtkMatrix4x4> matr;
	//d->InputVolumeNode->GetIJKToRASMatrix(matr);
	//d->OutputVolumeNode->SetIJKToRASMatrix(matr);

	double bounds29[6];
	inputVTK->GetBounds(bounds29);

	// Rotate about the center of the image
	vtkSmartPointer<vtkTransform> transform29 =
		vtkSmartPointer<vtkTransform>::New();

	// Compute the center of the image
	double center29[3];
	center29[0] = (bounds29[1] + bounds29[0]) / 2.0;
	center29[1] = (bounds29[3] + bounds29[2]) / 2.0;
	center29[2] = (bounds29[5] + bounds29[4]) / 2.0;

	// Rotate about the center
	transform29->Translate(center29[0], center29[1], center29[2]);
	transform29->RotateWXYZ((180+z), 0, 0, 1);
	transform29->Translate(-center29[0], -center29[1], -center29[2]);

	// Reslice does all of the work
	vtkSmartPointer<vtkImageReslice> reslice29 =
		vtkSmartPointer<vtkImageReslice>::New();
	reslice29->SetInputData(inputVTK);
	reslice29->SetResliceTransform(transform29);
	reslice29->SetInterpolationModeToCubic();
	reslice29->SetOutputSpacing(
		inputVTK->GetSpacing()[0],
		inputVTK->GetSpacing()[1],
		inputVTK->GetSpacing()[2]);
	reslice29->SetOutputOrigin(
		inputVTK->GetOrigin()[0],
		inputVTK->GetOrigin()[1],
		inputVTK->GetOrigin()[2]);
	reslice29->SetOutputExtent(inputVTK->GetExtent()); // Use a larger extent than the original image's to prevent clipping
	reslice29->Update();

	// Convert VTK -> ITK
	typedef itk::VTKImageToImageFilter< InputImageType > VTKtoITKFilterType;
	VTKtoITKFilterType::Pointer vtkToItkFilter = VTKtoITKFilterType::New();
	//vtkToItkFilter->SetInput(inputVTK);
	vtkToItkFilter->SetInput(reslice29->GetOutput());
	vtkToItkFilter->Update();
	InputImageType::Pointer image = vtkToItkFilter->GetOutput();

	// Software Guide : BeginLatex
	//
	// Creation of a \code{ResampleImageFilter} enables coordinates for
	// each of the pixels in the DRR image to be generated. These
	// coordinates are used by the \code{RayCastInterpolateImageFunction}
	// to determine the equation of each corresponding ray which is cast
	// through the input volume.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	using FilterType = itk::ResampleImageFilter<InputImageType, InputImageType >;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput(image);
	filter->SetDefaultPixelValue(0);
	//qDebug() << "function::projectVolume2()::ResampleImageFilter";

	// Software Guide : EndCodeSnippet
	// Software Guide : BeginLatex
	//
	// An Euler transformation is defined to position the input volume.
	// The \code{ResampleImageFilter} uses this transform to position the
	// output DRR image for the desired view.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	using TransformType = itk::CenteredEuler3DTransform< double >;
	TransformType::Pointer transform = TransformType::New();
	transform->SetComputeZYX(true);
	TransformType::OutputVectorType translation;
	translation[0] = 0;// tx; left -right
	translation[1] = 200;// ty; THIS WAS ORIGINALLY 0
	translation[2] = 0;// tz; up -down
	// constant for converting degrees into radians
	const double dtr = (std::atan(1.0) * 4.0) / 180.0;
	transform->SetTranslation(translation);
	//transform->SetRotation(dtr*rx, dtr*ry, dtr*rz);
	//transform->SetRotation(dtr*0., dtr*0., dtr*0.);
	transform->SetRotation(dtr*90., dtr*0., dtr*0.); //TODO change 3rd coordinate to 20
	//transform->SetRotation(dtr*(90.+x), dtr*y, dtr*z);
	//forward-back, clockwise-ccw, 


	InputImageType::PointType   imOrigin = image->GetOrigin();
	InputImageType::SpacingType imRes = image->GetSpacing();
	using InputImageRegionType = InputImageType::RegionType;
	using InputImageSizeType = InputImageRegionType::SizeType;
	InputImageRegionType imRegion = image->GetBufferedRegion();
	InputImageSizeType   imSize = imRegion.GetSize();
	imOrigin[0] += imRes[0] * static_cast<double>(imSize[0]) / 2.0;
	imOrigin[1] += imRes[1] * static_cast<double>(imSize[1]) / 2.0;
	imOrigin[2] += imRes[2] * static_cast<double>(imSize[2]) / 2.0;
	TransformType::InputPointType center;
	/*center[0] = cx + imOrigin[0];
	center[1] = cy + imOrigin[1];
	center[2] = cz + imOrigin[2];
	transform->SetCenter(center);*/
	center[0] = 0 + imOrigin[0];
	center[1] = 0 + imOrigin[1];
	center[2] = 0 + imOrigin[2];
	transform->SetCenter(center);
	// Software Guide : EndCodeSnippet
	// Software Guide : BeginLatex
	//
	// The \code{RayCastInterpolateImageFunction} is instantiated and passed the transform
	// object. The \code{RayCastInterpolateImageFunction} uses this
	// transform to reposition the x-ray source such that the DRR image
	// and x-ray source move as one around the input volume. This coupling
	// mimics the rigid geometry of the x-ray gantry.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	using InterpolatorType =
		itk::RayCastInterpolateImageFunction<InputImageType, double>;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	interpolator->SetTransform(transform);
	// Software Guide : EndCodeSnippet
	// Software Guide : BeginLatex
	//
	// We can then specify a threshold above which the volume's
	// intensities will be integrated.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet

	interpolator->SetThreshold(0);//threshold);

	// Software Guide : EndCodeSnippet
	// Software Guide : BeginLatex
	//
	// The ray-cast interpolator needs to know the initial position of the
	// ray source or focal point. In this example we place the input
	// volume at the origin and halfway between the ray source and the
	// screen. The distance between the ray source and the screen
	// is the "source to image distance" \code{sid} and is specified by
	// the user.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	InterpolatorType::InputPointType focalpoint;
	focalpoint[0] = imOrigin[0];
	focalpoint[1] = imOrigin[1];
	focalpoint[2] = imOrigin[2] - 400 / 2;//- sid / 2.;
	interpolator->SetFocalPoint(focalpoint);
	// Software Guide : EndCodeSnippet
	// Software Guide : BeginLatex
	//
	// Having initialised the interpolator we pass the object to the
	// resample filter.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	//interpolator->Print(std::cout);
	filter->SetInterpolator(interpolator);
	filter->SetTransform(transform);
	// Software Guide : EndCodeSnippet
	// Software Guide : BeginLatex
	//
	// The size and resolution of the output DRR image is specified via the
	// resample filter.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	// setup the scene
	InputImageType::SizeType   size;
	size[0] = 512;//dx;  // number of pixels along X of the 2D DRR image
	size[1] = 494;// dy;  // number of pixels along Y of the 2D DRR image
	size[2] = 1;   // only one slice
	filter->SetSize(size);
	InputImageType::SpacingType spacing;
	spacing[0] = 1;// sx;  // pixel spacing along X of the 2D DRR image [mm]
	spacing[1] = 1;// sy;  // pixel spacing along Y of the 2D DRR image [mm]
	spacing[2] = 1; // slice thickness of the 2D DRR image [mm]
	filter->SetOutputSpacing(spacing);
	// Software Guide : EndCodeSnippet

	// Software Guide : BeginLatex
	//
	// In addition the position of the DRR is specified. The default
	// position of the input volume, prior to its transformation is
	// half-way between the ray source and screen and unless specified
	// otherwise the normal from the "screen" to the ray source passes
	// directly through the centre of the DRR.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	double origin[dim];
	origin[0] = imOrigin[0] + 0 - 1.*((double)512 - 1.) / 2.;//o2Dx - sx*((double)dx - 1.) / 2.;
	origin[1] = imOrigin[1] + 0 - 1.*((double)494 - 1.) / 2.;//o2Dy - sy*((double)dy - 1.) / 2.;
	origin[2] = imOrigin[2] + 400 / 2.;//sid / 2.;
	filter->SetOutputOrigin(origin);
	// Software Guide : EndCodeSnippet
	// create writer
	//qDebug() << "function::projectVolume2()::RayCastImageFilter";

		// Software Guide : BeginLatex
		//
		// The output of the resample filter can then be passed to a writer to
		// save the DRR image to a file.
		//
		// Software Guide : EndLatex
		// Software Guide : BeginCodeSnippet
		using RescaleFilterType = itk::RescaleIntensityImageFilter<
			InputImageType, InputImageType >;
		RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
		rescaler->SetOutputMinimum(0);
		rescaler->SetOutputMaximum(255); //TODO USED TO BE 255
		rescaler->SetInput(filter->GetOutput());
		rescaler->Update();

		typedef itk::ImageFileWriter< InputImageType >  WriterType;
		WriterType::Pointer writer9 = WriterType::New();
		writer9->SetFileName("output_fixed.nrrd");
		writer9->SetInput(rescaler->GetOutput());
		writer9->Update();

   typedef itk::ExtractImageFilter< InputImageType, OutputImageType > SliceType3;
   SliceType3::Pointer slice3 = SliceType3::New();
   slice3->InPlaceOn();
   slice3->SetDirectionCollapseToSubmatrix();
   InputImageType::RegionType inputRegion3 =
 rescaler->GetOutput()->GetLargestPossibleRegion();
   InputImageType::SizeType size3 = inputRegion3.GetSize();
   size3[2] = 0;
   InputImageType::IndexType start3 = inputRegion3.GetIndex();
   const unsigned int sliceNumber3 = 0;
   start3[2] = sliceNumber3;
   InputImageType::RegionType desiredRegion3;
   desiredRegion3.SetSize(size3);
   desiredRegion3.SetIndex(start3);
   slice3->SetExtractionRegion(desiredRegion3);
   //qDebug() << "function::projectVolume2()::" << size3[0] << " " << size3[1] << " ~ " << start3[0] << " " << start3[1];
   slice3->SetInput(rescaler->GetOutput());
   slice3->Update();

		/*WriterType::Pointer writer4 = WriterType::New();
		writer4->SetFileName("output4.jpeg");

		using ExtractFilterType = itk::ExtractImageFilter< OutputImageType, OutputImageType >;
		ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
		extractFilter->SetDirectionCollapseToSubmatrix();

		// set up the extraction region [one slice]
		const OutputImageType * inputImage = rescaler->GetOutput();
		OutputImageType::RegionType inputRegion = inputImage->GetBufferedRegion();
		OutputImageType::SizeType size1 = inputRegion.GetSize();
		size1[2] = 1; // we extract along z direction
		OutputImageType::IndexType start = inputRegion.GetIndex();
		const unsigned int sliceNumber = 0;//atoi(argv[3]);
		start[2] = sliceNumber;
		OutputImageType::RegionType desiredRegion;
		desiredRegion.SetSize(size);
		desiredRegion.SetIndex(start);

		extractFilter->SetExtractionRegion(desiredRegion);
		using PasteFilterType = itk::PasteImageFilter< OutputImageType >;
		PasteFilterType::Pointer pasteFilter = PasteFilterType::New();
		using MedianFilterType = itk::MedianImageFilter< OutputImageType,
			OutputImageType >;
		MedianFilterType::Pointer medianFilter = MedianFilterType::New();
		extractFilter->SetInput(inputImage);
		medianFilter->SetInput(extractFilter->GetOutput());
		pasteFilter->SetSourceImage(medianFilter->GetOutput());
		pasteFilter->SetDestinationImage(inputImage);
		pasteFilter->SetDestinationIndex(start);

		OutputImageType::SizeType indexRadius;
		indexRadius[0] = 1; // radius along x
		indexRadius[1] = 1; // radius along y
		indexRadius[2] = 0; // radius along z
		medianFilter->SetRadius(indexRadius);
		medianFilter->UpdateLargestPossibleRegion();
		const OutputImageType * medianImage = medianFilter->GetOutput();
		pasteFilter->SetSourceRegion(medianImage->GetBufferedRegion());
		writer4->SetInput(pasteFilter->GetOutput());
		writer4->Update();*/

		/*OutputImageType::IndexType start1;
		start1[0] = 0;  // first index on X
		start1[1] = 0;  // first index on Y
		start1[2] = 0;  // first index on Z
		OutputImageType::SizeType  size1;
		size1[0] = 490;  // size along X
		size1[1] = 490;  // size along Y
		size1[2] = 1;  // size along Z
		OutputImageType::RegionType region;
		region.SetSize(size1);
		region.SetIndex(start1);

		typedef itk::ExtractImageFilter< OutputImageType, OutputImageType2 > ExtractImageFilterType;
		ExtractImageFilterType::Pointer extract = ExtractImageFilterType::New();
		extract->SetInput(rescaler->GetOutput());
		extract->SetDirectionCollapseToIdentity();
		extract->SetExtractionRegion(region);
		extract->Update();
		using WriterType4 = itk::ImageFileWriter< OutputImageType2 >;
		WriterType4::Pointer writer4 = WriterType4::New();
		writer4->SetFileName("output4.jpeg");
		writer4->SetInput(extract->GetOutput());*/

		using WriterType11 = itk::ImageFileWriter< OutputImageType >;
		WriterType11::Pointer writer11 = WriterType11::New();
		writer11->SetFileName("output_fixed2.nrrd");
		writer11->SetInput(slice3->GetOutput());
		writer11->Update();


		using ITKtoVTKFilterType = itk::ImageToVTKImageFilter< OutputImageType >;
		ITKtoVTKFilterType::Pointer itktovtkfilter = ITKtoVTKFilterType::New();
		itktovtkfilter->SetInput(slice3->GetOutput());
		itktovtkfilter->Update();
		vtkImageData * outputVTK = itktovtkfilter->GetOutput();

		/*vtkSmartPointer<vtkImageReslice> reslice =
			vtkSmartPointer<vtkImageReslice>::New();
		reslice->SetOutputExtent(0, 9, 0, 100, 0, 0);
		reslice->SetInputData(outputVTK);
		reslice->Update();*/

		/*vtkSmartPointer<vtkImageWriter> writer5 =
			vtkSmartPointer<vtkImageWriter>::New();
		writer5->SetInputData(outputVTK);
		writer5->SetFileName("output3.vtk");
		writer5->Write();*/

		vtkSmartPointer<vtkImageCast> castFilter =
			vtkSmartPointer<vtkImageCast>::New();
		castFilter->SetOutputScalarTypeToUnsignedChar();
		castFilter->SetInputData(outputVTK);
		castFilter->Update();

		d->OutputVolumeNode->SetAndObserveImageData(castFilter->GetOutput());

		const auto duration = std::chrono::duration_cast<ms>(clock::now() - before);

		std::cout << "It took " << duration.count() / 1000.0 << "s" << std::endl;
}

void qSlicerVertRegModuleWidget::cropVolume()
{
	Q_D(qSlicerVertRegModuleWidget);

	using clock = std::chrono::system_clock;
	using ms = std::chrono::milliseconds;

	const auto before = clock::now();

	/*vtkSlicerApplicationLogic* slicerAppLogic = this->module()->appLogic();
	vtkMRMLSelectionNode* selectionNode = slicerAppLogic->GetSelectionNode();
	selectionNode->SetReferenceActiveVolumeID(volumeNode->GetID());
	slicerAppLogic->PropagateVolumeSelection();*/

	qDebug() << "function::segmentVolume()";

	const int dim = 3;

	typedef short PType;
	typedef itk::Image< PType, dim > InputImageType;
	typedef itk::Image< PType, dim > OutputImageType;


	vtkSmartPointer<vtkImageData> inputVTK = vtkSmartPointer<vtkImageData>::New();
	inputVTK = d->InputVolumeNode->GetImageData();

	//vtkNew<vtkMatrix4x4> matr;
	//d->InputVolumeNode->GetIJKToRASMatrix(matr);
	//d->OutputVolumeNode->SetIJKToRASMatrix(matr);

	qDebug() << "function::projectVolumeFilter1()";

	typedef itk::VTKImageToImageFilter< InputImageType > VTKtoITKFilterType;
	VTKtoITKFilterType::Pointer vtkToItkFilter = VTKtoITKFilterType::New();
	vtkToItkFilter->SetInput(inputVTK);
	vtkToItkFilter->Update();
	InputImageType::Pointer image = vtkToItkFilter->GetOutput();

	qDebug() << "function::projectVolumeFilter2()";


	/*typedef itk::MeanProjectionImageFilter< InputImageType, OutputImageType > ProjectionFilterType;
	ProjectionFilterType::Pointer projfilter = ProjectionFilterType::New();
	projfilter->SetInput(image);
	projfilter->SetProjectionDimension(1);
	projfilter->Update();*/

	using FilterType = itk::ThresholdImageFilter< InputImageType >;
	FilterType::Pointer thresholdfilter = FilterType::New();
	thresholdfilter->SetInput(image);
	thresholdfilter->ThresholdOutside(95, 500);
	thresholdfilter->SetOutsideValue(0);
	thresholdfilter->Update();

	//itk::SimpleFilterWatcher watcher(projfilter, "filter");
	qDebug() << "function::projectVolumeFilter3()";


	OutputImageType::Pointer outputITK = thresholdfilter->GetOutput();

	qDebug() << "function::projectVolumeFilter4()";

	/*typedef itk::VTKImageExport<InputImageType> ImageExportType;

	// Create the itk::VTKImageExport instance and connect it to the
	// itk::CurvatureFlowImageFilter.
	ImageExportType::Pointer itkExporter = ImageExportType::New();
	itkExporter->SetInput(outputITK);

	// Create the vtkImageImport and connect it to the
	// itk::VTKImageExport instance.
	vtkImageImport* vtkImporter = vtkImageImport::New();
	ConnectPipelines(itkExporter, vtkImporter);
	vtkImporter->Update();

	vtkImageShiftScale* shifter = vtkImageShiftScale::New();

	shifter->SetInputData(vtkImporter->GetOutput());
	shifter->SetScale(256);
	shifter->SetOutputScalarTypeToUnsignedChar();*/



	qDebug() << "function::projectVolumeFilter5()";

	/*using ReaderType = itk::ImageFileReader< OutputImageType >;
	ReaderType::Pointer reader2 = ReaderType::New();
	reader2->SetFileName("output_threshold.nrrd");
	reader2->Update();*/

	using ITKtoVTKFilterType2 = itk::ImageToVTKImageFilter< OutputImageType >;
	ITKtoVTKFilterType2::Pointer itktovtkfilter2 = ITKtoVTKFilterType2::New();
	itktovtkfilter2->SetInput(outputITK);
	try
	{
		itktovtkfilter2->Update();
	}
	catch (itk::ExceptionObject & error)
	{
		std::cerr << "Error: " << error << std::endl;
	}
	//vtkSmartPointer<vtkImageData> outputVTK = vtkSmartPointer<vtkImageData>::New();
	//outputVTK = itktovtkfilter2->GetOutput();
	vtkImageData * outputVTK = itktovtkfilter2->GetOutput();


	//outputVTK->Print(std::cout);
	//vtkImageData * outputVTK = itktovtkfilter->GetOutput();

	//outputVTK = itktovtkfilter2->GetOutput();

	qDebug() << "function::projectVolumeFilter6()";


		/*vtkSmartPointer<vtkImageReader> reader = vtkSmartPointer<vtkImageReader>::New();
		//writer->SetInputData(outputVTK);
		reader->SetInputConnection()
		reader->Update();
		outputVTK = reader->GetOutput();*/

		

	/*vtkSmartPointer<vtkImageCast> castFilter2 =
	vtkSmartPointer<vtkImageCast>::New();
	castFilter2->SetOutputScalarTypeToShort();
	castFilter2->SetInputData(outputVTK);
	castFilter2->Update();*/

	//d->OutputVolumeNode->SetAndObserveImageData(castFilter->GetOutput());
	d->OutputVolumeNode->SetAndObserveImageData(outputVTK);

	const auto duration = std::chrono::duration_cast<ms>(clock::now() - before);

	std::cout << "It took " << duration.count() / 1000.0 << "ms" << std::endl;
}

void qSlicerVertRegModuleWidget::rotateVolume()
{
	Q_D(qSlicerVertRegModuleWidget);

	using clock = std::chrono::system_clock;
	using ms = std::chrono::milliseconds;

	const auto before = clock::now();

	qDebug() << "function::rotateVolume()";

	const int dim = 3;

	typedef short PType;
	typedef itk::Image< PType, dim > InputImageType;
	typedef itk::Image< PType, dim > OutputImageType;


	vtkSmartPointer<vtkImageData> inputVTK = vtkSmartPointer<vtkImageData>::New();
	inputVTK = d->InputVolumeNode->GetImageData();

	vtkNew<vtkMatrix4x4> matr;
	d->InputVolumeNode->GetIJKToRASMatrix(matr);
	d->OutputVolumeNode->SetIJKToRASMatrix(matr);

	/*vtkSmartPointer<vtkImageCast> castFilter2 =
		vtkSmartPointer<vtkImageCast>::New();
	castFilter2->SetOutputScalarTypeToUnsignedChar();
	castFilter2->SetInputData(inputVTK);
	castFilter2->Update();*/

		typedef itk::VTKImageToImageFilter< InputImageType > VTKtoITKFilterType;
		VTKtoITKFilterType::Pointer vtkToItkFilter = VTKtoITKFilterType::New();
		vtkToItkFilter->SetInput(inputVTK);
		vtkToItkFilter->Update();
		InputImageType::Pointer image = vtkToItkFilter->GetOutput();
	
		typedef itk::CenteredAffineTransform< double , dim > TransformType;
		typedef itk::ResampleImageFilter<InputImageType, InputImageType> RotateFilterType;
		RotateFilterType::Pointer rotatefilter = RotateFilterType::New();
		TransformType::Pointer transform = TransformType::New();

		TransformType::InputPointType origin;
		InputImageType::SizeType size;


		rotatefilter->SetTransform(transform);
		size = image->GetLargestPossibleRegion().GetSize();

		rotatefilter->SetInput(image);
		TransformType::InputPointType center;

		origin = image->GetOrigin();
		center[0] = origin[0] + size[0] / 2.0;
		center[1] = origin[1] + size[1] / 2.0;
		center[2] = origin[2] + size[2] / 2.0;

		transform->SetCenter(center);
		const double degreesToRadians = atan(1.0) / 45.0;
		const double angle = 30 * degreesToRadians;
		TransformType::OutputVectorType axis;
		axis[0] = 0;
		axis[1] = 0;
		axis[2] = 1;

		//transform->ComputeOffset();
		transform->Rotate3D(axis, angle, false);

		rotatefilter->SetDefaultPixelValue(100);
		rotatefilter->SetSize(size);
		rotatefilter->Update();

		if (debug) {
			/*typedef itk::Image< unsigned char, dim - 1 > WriterImageType;
			using CastFilterType = itk::CastImageFilter<OutputImageType, WriterImageType >;
			CastFilterType::Pointer castFilteritk = CastFilterType::New();
			castFilteritk->SetInput(outputITK);
			castFilteritk->Update();

			typedef itk::ImageFileWriter< WriterImageType > WriterType;*/
			typedef itk::ImageFileWriter< OutputImageType > WriterType;
			WriterType::Pointer writer = WriterType::New();
			//writer->SetInput(castFilteritk->GetOutput());
			writer->SetInput(rotatefilter->GetOutput());
			writer->SetFileName("output_rotate.nrrd");
			writer->Update();
		}

		/*using ITKtoVTKFilterType = itk::ImageToVTKImageFilter< OutputImageType >;
		ITKtoVTKFilterType::Pointer itktovtkfilter = ITKtoVTKFilterType::New();
		itktovtkfilter->SetInput(rotatefilter->GetOutput());
		itktovtkfilter->Update();*/
		//vtkImageData * outputVTK = itktovtkfilter->GetOutput();
		//outputVTK->Print(std::cout);
		vtkSmartPointer<vtkImageData> outputVTK = vtkSmartPointer<vtkImageData>::New();

		if (debug) {

			vtkSmartPointer<vtkNrrdReader> writer = vtkSmartPointer<vtkNrrdReader>::New();
			writer->SetInputData(outputVTK);
			writer->SetFileName("output_rotate.nrrd");
			writer->Update();
			outputVTK = writer->GetOutput();
		}

		/*vtkSmartPointer<vtkImageReslice> reslice =
			vtkSmartPointer<vtkImageReslice>::New();
		reslice->SetOutputExtent(0, 9, 0, 100, 0, 0);
		reslice->SetInputData(outputVTK);
		reslice->Update();*/

		/*vtkSmartPointer<vtkImageWriter> writer5 =
			vtkSmartPointer<vtkImageWriter>::New();
		writer5->SetInputData(outputVTK);
		writer5->SetFileName("output3.vtk");
		writer5->Write();*/

		/*vtkSmartPointer<vtkImageCast> castFilter =
			vtkSmartPointer<vtkImageCast>::New();
		castFilter->SetOutputScalarTypeToShort();
		castFilter->SetInputData(outputVTK);
		castFilter->Update();*/

		/*vtkSmartPointer<vtkPNGWriter> writer6 =
			vtkSmartPointer<vtkPNGWriter>::New();
		writer6->SetFileName("output6.png");
		writer6->SetInputData(castFilter->GetOutput());
		writer6->Write();*/

		d->OutputVolumeNode->SetAndObserveImageData(outputVTK);

		const auto duration = std::chrono::duration_cast<ms>(clock::now() - before);

		std::cout << "It took " << duration.count() / 1000.0 << "ms" << std::endl;
}


//-----------------------------------------------------------------------------
//-------------Ignore this method----------------------------------------------
void qSlicerVertRegModuleWidget::projectVolume()
{
	Q_D(qSlicerVertRegModuleWidget);

	qDebug() << "function::projectVolume()";
	//qDebug() << d->InputVolumeNode->GetID();

	const int dim = 3;

	typedef short PType;
	typedef itk::Image< PType, dim > ImageType;
	typedef itk::Image< unsigned char, 3 > OutputImageType;

	vtkSmartPointer<vtkImageData> inputVTK = vtkSmartPointer<vtkImageData>::New();
	inputVTK = d->InputVolumeNode->GetImageData();

	
	vtkNew<vtkMatrix4x4> matr;
	d->InputVolumeNode->GetIJKToRASMatrix(matr);
	d->OutputVolumeNode->SetIJKToRASMatrix(matr);
	//d->OutputVolumeNode->SetAndObserveImageData(inputVTK);
	

	/*vtkNew<vtkMatrix4x4> matr;
	d->InputVolumeNode->GetIJKToRASMatrix(matr);
	d->InputVolumeNode->SetRASToIJKMatrix(matr);*/

	/*vtkNew<vtkMatrix4x4> matr;
	vtkNew<vtkMatrix4x4> matrinv;
	d->InputVolumeNode->GetRASToIJKMatrix(matr);
	vtkMatrix4x4::Invert(matr, matrinv);
	d->InputVolumeNode->SetRASToIJKMatrix(matrinv);*/

	

	typedef itk::VTKImageToImageFilter< ImageType > VTKtoITKFilterType;
	VTKtoITKFilterType::Pointer vtkToItkFilter = VTKtoITKFilterType::New();
	vtkToItkFilter->SetInput(inputVTK);
	vtkToItkFilter->Update();
	ImageType::Pointer inputITK = vtkToItkFilter->GetOutput();


	//typedef itk::ImageFileWriter< ImageType > WriterTyp4e;
	//WriterTyp4e::Pointer writer4 = WriterTyp4e::New();
	//writer4->SetInput(inputITK);
	//writer4->SetFileName("output3.dcm");
	//writer4->Update();

	/*itk::FixedArray<bool, 3> flipAxes;
	flipAxes[0] = false;
	flipAxes[1] = true;
	flipAxes[2] = false;
	typedef itk::FlipImageFilter <ImageType> FlipImageFilterType;
	FlipImageFilterType::Pointer flipFilter = FlipImageFilterType::New();
	flipFilter->SetInput(inputITK);
	flipFilter->SetFlipAxes(flipAxes);
	ImageType::Pointer inputITK2 = flipFilter->GetOutput();*/

	//vtkMatrix4x4 *mat = new vtkMatrix4x4();
	//vtkSmartPointer<vtkMatrix4x4> mat = vtkSmartPointer<vtkMatrix4x4>::New();




	// Convert from LPS (ITK) to RAS (Slicer)
	// input: transformVtk_LPS matrix in vtkMatrix4x4 in resampling convention in LPS
	// output: transformVtk_RAS matrix in vtkMatri4x4 in modeling convention in RAS

	/*double * origin = d->InputVolumeNode->GetOrigin();
	origin[0] = -origin[0];
	origin[1] = -origin[1];
	origin[2] = -origin[2];
	double * spacing = d->InputVolumeNode->GetSpacing();
	inputITK->SetOrigin(origin);
	inputITK->SetSpacing(spacing);*/
	//inputITK->SetSpacing(d->InputVolumeNode->GetSpacing());



	typedef itk::MeanProjectionImageFilter< ImageType, OutputImageType > ProjectionFilterType;
	ProjectionFilterType::Pointer projfilter = ProjectionFilterType::New();
	projfilter->SetInput(inputITK);
	projfilter->SetProjectionDimension(1);

	itk::SimpleFilterWatcher watcher(projfilter, "filter");

	typedef itk::ImageFileWriter< OutputImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	OutputImageType::Pointer outputITK = projfilter->GetOutput();

	typedef itk::ImageFileWriter< OutputImageType > WriterTyp3e;
	WriterTyp3e::Pointer writer3 = WriterTyp3e::New();
	writer3->SetInput(projfilter->GetOutput());
	writer3->SetFileName("output2.dcm");
	writer3->Update();

	using ITKtoVTKFilterType = itk::ImageToVTKImageFilter< OutputImageType >;
	ITKtoVTKFilterType::Pointer itktovtkfilter = ITKtoVTKFilterType::New();
	itktovtkfilter->SetInput(projfilter->GetOutput());
	itktovtkfilter->Update();
	//vtkImageData * outputVTK = itktovtkfilter->GetOutput();
	vtkSmartPointer<vtkJPEGReader> vtkimagefilter =
		vtkSmartPointer<vtkJPEGReader>::New();
	vtkimagefilter->SetFileName("output2.jpeg");
	vtkimagefilter->Update();
	vtkImageData * outputVTK = vtkimagefilter->GetOutput();



	d->OutputVolumeNode->SetAndObserveImageData(outputVTK);

}




//typedef itk::ImageFileWriter< ImageType > WriterType;
//WriterType::Pointer writer = WriterType::New();
//ImageType::Pointer outputITK = projfilter->GetOutput();
//using ITKtoVTKFilterType = itk::ImageToVTKImageFilter< ImageType >;
//ITKtoVTKFilterType::Pointer itktovtkfilter = ITKtoVTKFilterType::New();
//itktovtkfilter->SetInput(outputITK);
//itktovtkfilter->Update();
//vtkImageData * outputVTK = itktovtkfilter->GetOutput();
////d->OutputVolumeNode->SetOrigin(d->InputVolumeNode->GetOrigin());
//d->OutputVolumeNode->SetAndObserveImageData(outputVTK);


//writer->SetImageIO(d->OutputVolumeNode);
//MRMLIDImageIO* example = new MRMLIDImageIO();
//itk::MRMLIDImageIO* example = new MRMLIDImageIO();
//typedef itk::MRMLIDImageIO mrmlIO;
//mrmlIO::Pointer outputmrml = mrmlIO::New();
//vtkMRMLNode* node = (vtkMRMLNode*)d->OutputVolumeNode;

//vtkMRMLStorableNode* storable = vtkMRMLStorableNode::SafeDownCast(node);
//vtkMRMLNode* node = (vtkMRMLNode*)d->OutputVolumeNode;

//vtkMRMLStorableNode* storable = vtkMRMLStorableNode::SafeDownCast(node);
//itk::MRMLIDImageIO * f = new itk::MRMLIDImageIO::New();
//writer->SetFileName(d->OutputVolumeNode.c_str);
//writer->SetImageIO(itk::MRMLIDImageIO::New());
//writer->SetImageIO(f);
//writer->SetFileName("output.mrml");
//writer->SetFileName("testfilenameplz");
//writer->SetImageIO(outputmrml);
//writer->Update();
//
//void qSlicerVertRegModuleWidget::setMRMLScene(vtkMRMLScene* scene)
//{
//	this->Superclass::setMRMLScene(scene);
//	
//	// observe close event so can re-add a parameters node if necessary
//	/*qvtkReconnect(this->mrmlScene(), vtkMRMLScene::EndImportEvent, this, SLOT(onMRMLSceneEndBatchProcessEvent()));
//	qvtkReconnect(this->mrmlScene(), vtkMRMLScene::EndBatchProcessEvent, this, SLOT(onMRMLSceneEndBatchProcessEvent()));
//	qvtkReconnect(this->mrmlScene(), vtkMRMLScene::EndCloseEvent, this, SLOT(onMRMLSceneEndBatchProcessEvent()));
//	qvtkReconnect(this->mrmlScene(), vtkMRMLScene::EndRestoreEvent, this, SLOT(onMRMLSceneEndBatchProcessEvent()));
//	this->updateWidgetFromMRML();*/

//qSlicerCoreIO::

//reader->SetFileName(d->InputVolumeNode->GetID());
/*itk::ImageIOBase::Pointer imageIO =
itk::ImageIOFactory::CreateImageIO(
d->InputVolumeNode->GetID, itk::ImageIOFactory::ReadMode);*/
//itk::MRMLIDImageIO::Pointer imageIO = itk::MRMLIDImageIOFactory::FactoryNew();
//itk::ObjectFactoryBase::RegisterFactory(itk::MRMLIDImageIOFactory::New());
// Read the file
//reader2 = vtkSmartPointer< vtkITKArchetypeImageSeriesScalarReader >::New();
//reader2->SetArchetype(d->InputVolumeNode->GetID());
//reader2->SetOutputScalarTypeToNative();
//reader2->SetDesiredCoordinateOrientationToNative();
//reader2->SetUseNativeOriginOn();
//reader2->Update();
//vtkSmartPointer< vtkImageChangeInformation > ici = vtkSmartPointer< vtkImageChangeInformation >::New();
//ici->SetInput(reader->GetOutput());
//ici->SetOutputSpacing(1, 1, 1);
//ici->SetOutputOrigin(0, 0, 0);
//ici->Update();

//image = ici->GetOutput();
//image->Update();
//d->InputVolumeNode->GetImageData();
//ImageIOType::Pointer vtkIO = ImageIOType::New(); d->InputVolumeNode->GetImageData())
//typedef itk::MRMLIDImageIO ImageIOType;

//using ImageIOType = itk::MRMLIDImageIO;
////typedef itk::ImageFileReader< ImageIOType > ReaderType2;
//typedef itk::Image< PType, dim > InputImageType;
//using ReaderType2 = itk::ImageFileReader< InputImageType >;
//ReaderType2::Pointer reader2 = ReaderType2::New();

//reader2->SetFileName(d->InputVolumeNode->GetID());
//reader2->Update();
//}
//itk::VTKImageIO::Pointer vtkImgIo = itk::VTKImageIO::New();
//this->mrmlScene;
//
//vtkImgIo->SetFileName("slicer:<scene id>#<node id>");
//vtkImgIo->Read(imageITK);

//reader->SetImageIO(vtkImgIo);

//ImageIOType::Pointer gdcmIO = new ImageIOType();
//itk::ObjectFactoryBase::RegisterFactory(itk::MRMLIDImageIOFactory::New());
//ImageIOType mrmlo = itk::ObjectFactoryBase::CreateInstance("MRMLIDImageIO");
//ImageIOType::Pointer imIO =
//itk::MRMLIDImageIOFactory::CreateInstance(itk::MRMLIDImageIO);
//itk::MRMLFactory::CreateImageIO(
//	inputFilename.c_str(), itk::ImageIOFactory::ReadMode);
//reader2->SetImageIO(itk::MRMLIDImageIO::New());
//itk::LightObject::Pointer mrmlop = itk::ObjectFactoryBase::CreateInstance("MRMLIDImageIO");
//qDebug() << mrmlop.Print();
//reader2->GetOutput();
//reader2->SetImageIO(mrmlop);
//itk::MRMLIDImageIOFactory::Read

//reader->SetImageIO(imageIO);
//reader->SetFileName(d->InputVolumeNode->GetID());
//reader->SetFileName("D:\Med school\IR Research\QIN-HEADNECK\QIN-HEADNECK-01-0021\10-20-1987-PET1PETCTHeadNeck Adult-63581\20061221-CT HeadNeck  3.0  B30fABDOMEN-24378/000000.dcm");
//reader->SetImageIO(gdcmIO);
//------------------------------------------------------------
//typedef itk::MeanProjectionImageFilter< IType, IType > FilterType;
//FilterType::Pointer filter = FilterType::New();
//filter->SetInput(reader->GetOutput());

//itk::SimpleFilterWatcher watcher(filter, "filter");

//typedef itk::ImageFileWriter< IType > WriterType;
//WriterType::Pointer writer = WriterType::New();
//writer->SetInput(filter->GetOutput());
//------------------------------------------------------------

//typedef unsigned char PType;
//typedef itk::Image< PType, dim > IType;

//typedef itk::ImageFileReader< IType > ReaderType;
//ReaderType::Pointer reader = ReaderType::New();
