from pydantic import BaseModel, PositiveFloat, Field, model_validator
from typing import List, Optional, Tuple, Annotated, ClassVar
from typing_extensions import Self
import numpy as np

EffMax = Annotated[float, Field(ge=0.0, le=1.0)]
EffSig = Annotated[float, Field(ge=0.0, le=4.0)]
EffM0 = Annotated[float, Field(ge=20.0, le=31.0)]


class CharacterizationParameter(BaseModel, validate_assignment=True):
    """
    Base class for all efficiency file parameters
    """
    name: ClassVar[str] = 'base'

    def __str__(self):
        values = " ".join([ f"{self.__dict__.get(key, None)}"
                  for key in self.eff_file_format_order ])
        return f"{self.name}= {values}"

    @classmethod
    def from_string(cls, line: str) -> Self:
        """
        Create a parameter object from a string
        """
        keys = cls.eff_file_format_order
        values = line.split('=')[1].strip().split()
        return cls(**{key: value for key, value in zip(keys, values)})

    @property
    def eff_file_format_order(self) -> List[str]:
        """
        Return list of parameter names in order needed for .eff file format
        """
        return self.__pydantic_fields__.keys()


class Bandpass(CharacterizationParameter):
    name: ClassVar[str] = 'filter'
    Bandpass: str = Field(default='r', pattern='[ugriz]')


class RateCut(CharacterizationParameter):
    # The rate cut is defined by a minimum and maximum rate
    # The rate is defined in arcseconds per hour
    name: ClassVar[str] = 'rate_cut'
    rate_min: PositiveFloat = Field(ge=0.0)
    rate_max: PositiveFloat = Field(ge=0.0)

    @model_validator(mode='after')
    def check_rate(self) -> Self:
        """
        Check that rate_min is less than rate_max
        """
        if self.rate_min >= self.rate_max:
            raise ValueError('rate_max must be greater than rate_min')
        return self


class PhotFraction(CharacterizationParameter):
    # The fraction of objects with 1, 2 or 3 photometric measurements
    # use to determine the official magnitude
    name : ClassVar[str] = 'phot_frac'
    frac_one: float = Field(ge=0.0, le=1.0)
    frac_two: float = Field(ge=0.0, le=1.0)
    frac_three: float = Field(ge=0.0, le=1.0)

    @model_validator(mode='after')
    def check_phot_frac(self) -> Self:
        """
        Check that the sum of the photometric fractions is 1
        """
        if self.frac_one + self.frac_two + self.frac_three != 1.0:
            raise ValueError('The sum of the photometric fractions must be 1')
        return self
    

class TrackingFraction(CharacterizationParameter):
    name: ClassVar[str] = 'track_frac'
    peak: float = Field(1.0, ge=0.0, le=1.0)
    Rc: float = Field(24.0, ge=22.0, le=31.0)
    slope: float = Field(-5, ge=-10.0, le=0.0)  # starting at Rc drop from peak to 0 with slope


class MagError(CharacterizationParameter):
    # The magnitude error is defined by a set of parameters
    # that are used to define the growth of the uncertainty
    name: ClassVar[str] = 'mag_error'
    magerr_bright: float = 0.01  # magerr for mag < 21
    magerr_slope: float = 0.2 # magerr slope for 21 < mag < mag_mid
    mag_mid: float = 24.0 # range until which the magerr slope is defined
    magerr_faint_slope: float = 0.5 # slope on grown of magerr for mag > mag_mid
    mag_faint: float = 24.2 # mag to start the faint magerr slope
    magerr_bias: float = -0.2 # bias for mag > mag_faint

    def __call__(self, mag) -> np.ndarray:
        """
        Simulate the magnitude error as a function of input magnitude.
        This method calculates the magnitude error (`merr`) based on the input magnitude (`mag`) 
        and the object's internal parameters. It applies different error models depending on 
        the brightness of the input magnitude and introduces random scatter and bias adjustments.
        Parameters:
        -----------
        mag : np.ndarray
            The input magnitude(s) for which the uncertainty estimate is needed.
        Returns:
        --------
        np.ndarray
            The simulated magnitude error for the input magnitudes.
        """

        mag_th = mag
        magerr_bright = self.magerr_bright
        magerr_mid =    self.magerr_bright*10.0**(self.magerr_slope*(mag_th   -21.0))
        magerr_faint =  self.magerr_bright*10.0**(self.magerr_slope*(self.mag_mid-21.0)) - (mag_th - self.mag_mid)*self.magerr_faint_slope
        magerr_bias = self.magerr_bias
        magerr_faint[magerr_faint < 0] = 0

        tmp = np.random.random_sample(len(mag_th))
        A = np.sqrt(6.0)*(np.sqrt(2*tmp) - 1)
        B = np.sqrt(6.0)*(1 - np.sqrt(2*(1-tmp)))
        tmp[tmp <= 0.5] = A[tmp <= 0.5]
        tmp[tmp > 0.5] = B[tmp > 0.5]
        
        magerr = magerr_bright*(mag <= 21.0) + magerr_mid*((mag > 21.0) & ( mag <= self.mag_mid)) + magerr_faint*(mag > self.mag_mid)
        # the bias slope is defined for the mag after the random scatter is defined.
        _mag = mag_th + magerr*tmp 
        _mag[mag_th > self.mag_faint] += (_mag[mag_th > magerr_faint] - self.mag_faint)*self.magerr_bias
        return mag-_mag


class LinearParam(CharacterizationParameter):
    """
    LinearParam class to define the parameters for the linear function

    :math:`efficiency(m) = A * (m - R_1) / (R_2 - R_1)`
    where :math:`A` is the efficiency at :math:`m = R_1`, and :math:`R_1` and
    :math:`R_2` are the two points where the efficiency is defined. The
    efficiency is 0 at :math:`m = R_2` and :math:`A` at :math:`m = R_1`.
    The efficiency is defined as:
    .. math::
        \begin{eqnarray}
        A & {\rm if} & m < R_1 \\
        \frac{(m - R_1) A}{R_2 - R_1} & {\rm if} & R_1 \le m < R_2 \\
        0 & {\rm if} & m \ge R_2
    \end{eqnarray}
    """
    # The parameters are defined as a list of floats
    name: ClassVar[str] = 'linear_param'
    A: EffMax = 1.0
    R_1: EffM0 = 24.0
    R_2: EffM0 = 26.0


class DoubleTanhParam(CharacterizationParameter):
    """
    DoubleTanh class to define the parameters for the double tanh function

    :math:`efficiency(m) = A * (tanh((m - R_1)/w_1) + tanh((R_2 - m)/w_2))`
    where :math:`A` is the efficiency at :math:`m = R_1`, and :math:`R_1` and
    :math:`R_2` are the two points where the efficiency is defined. The
    efficiency is 0 at :math:`m = R_2` and :math:`A` at :math:`m = R_1`.

    """
    # The parameters are defined as a list of floats
    name: ClassVar[str] = 'double_param'
    A: EffMax = 1.0
    M_0: EffM0 = 24.0
    w_1: EffSig = 1.0
    w_2: EffSig = 0.5

    
class SquareParam(CharacterizationParameter):
    """
    SquareEff class to define the parameters for the square function

    (eff_max - c * (R-21.0)^2) / (1. + exph((R-M_0)/sig))
    """
    # The parameters are defined as a list of floats
    name: ClassVar[str] = 'square_param'
    eff_max: EffMax =  1.0
    c: EffSig = 0.05
    M_0: EffM0 = 24.0
    sig: EffSig = 0.5


class LookupParam(CharacterizationParameter):
    """
    LookupEff class to define the parameters for the lookup function
    """
    name: ClassVar[str] = 'lookup_param'
    # The parameters are defined as a list of floats
    mag: float
    eff: float


class Function(CharacterizationParameter):
    name: ClassVar[str] = 'function'
    function: str = Field(default='double')


class Rates(CharacterizationParameter):
    """
    RateEfficiency class to define the detection efficiency at a particular magnitude given the rate of motion
    """
    # The efficiency is defined over rate ranges
    name: ClassVar[str] = 'rates'
    rate_min: PositiveFloat
    rate_max: PositiveFloat

    @model_validator(mode='after')
    def check_rate(self) -> Self:
        """
        Check that rate_min is less than rate_max
        """
        if self.rate_min >= self.rate_max:
            raise ValueError('rate_max must be greater than rate_min')
        return self


class MagLim(CharacterizationParameter):
    """
    MagLim class to define the limiting magnitude of the block
    """
    name: ClassVar[str] = 'mag_lim'
    mag_lim: EffM0 = Field(default=24.0, ge=20.0, le=31.0)
    # The limiting magnitude is defined as a float


class ParamFactory(object):
    """
    Class to define the param record for efficiency function
    """

    def __init__(self, base_class=CharacterizationParameter):
        self._registry = { }
        _ = [ self.register(cls.name, cls) for cls in base_class.__subclasses__() ]

    def register(self, name: str, cls: type):
        self._registry[name] = cls

    def __call__(self, name: str):
        return self._registry[name]


class Characterization(BaseModel):
    """
    SearchRates class to define the search rates for a given efficiency file
    """
    components: List[CharacterizationParameter]

    @classmethod
    def from_eff_file(cls, filename: str) -> Self:
        """
        Create a SearchRates object from a string
        """
        with open(filename, 'r') as fojb:
            lines = fojb.readlines()
        # Remove comments and empty lines   
        lines = [line for line in lines if not line.startswith('#') and line.strip() != '']
        components = []
        for line in lines:
            component = line.split('=')[0].strip()
            components.append(ParamFactory(component).from_string(line))
        return cls(components=components)

    def __str__(self):
        """
        String representation of the SearchRates class
        """
        value = ""
        for component in self.components:
            value += f"{component}\n"
        return value

function = 'double'
M_0 = 24.0
A = 0.9
w1 = 0.5
w2 = 0.5

characterization = Characterization(
    components=[
    RateCut(rate_min=0.5, rate_max=10.0),
    MagError(magerr_bright=0.2,
                       magerr_slope=0.5,
                       mag_mid=M_0+w1,
                       magerr_faint_slope=0.5,
                       mag_faint=M_0+w2,
                       magerr_bias=-0.2),
    PhotFraction(frac_one=0.5, frac_two=0.3, frac_three=0.2),
    TrackingFraction(
        peak=1.0,
        Rc=M_0+w2,
        slope=-5.0
    ),
    Bandpass(Bandpass='r'),
    Rates(rate_min=0.5, rate_max=10.0),
    Function(function='double'),
    DoubleTanhParam(A=1.0, M_0=24.0, w_1=1.0, w_2=0.5),
    MagLim(mag_lim=24.0)
    ])

print(characterization)
