
import { Route, Routes } from 'react-router-dom';


  const AppRoutes = ({ pages }) => (
    <Routes>
      <Route path='/' element={pages['Home']}></Route>
      {Object.keys(pages).map((page) => (
        <Route path={`/${page.toLowerCase()}`} element={pages[page]} />
      ))}
    </Routes>
  );

  export default AppRoutes;